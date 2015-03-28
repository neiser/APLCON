#include "APLCON.hpp"

// detail code is in namespace APLCON_ (note the underscore)
#include "detail/APLCON_cc.hpp"

// long ostream stuff is in extra header
#include <detail/APLCON_ostream.hpp>

// include the wrapping C functions
extern "C" {
#include "wrapper/APLCON.h"
}

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>


using namespace std;

int APLCON::instance_counter = 0;
int APLCON::instance_lastfit = 0;

// short NaN is really handy...
const double APLCON::NaN = numeric_limits<double>::quiet_NaN();

const APLCON::Variable_Settings_t APLCON::Variable_Settings_t::Default = {
  APLCON::Distribution_t::Gaussian,
  {
    -numeric_limits<double>::infinity(),
    numeric_limits<double>::infinity()
  },
  APLCON::NaN
};

// transferred to APLCON in Init() method
const APLCON::Fit_Settings_t APLCON::Fit_Settings_t::Default = {
  0,  // no debug printout from APLCON
  -1, // max iterations
  APLCON::NaN,
  APLCON::NaN,
  APLCON::NaN,
  APLCON::NaN
};

// proper default result
const APLCON::Result_t APLCON::Result_t::Default = {
  "", // Name
  APLCON::Result_Status_t::_Unknown,
  APLCON::NaN,
  -1,
  APLCON::NaN,
  -1,
  -1,
  {},
  {}
};

// method implementations

void APLCON::AddMeasuredVariable(const std::string &name, const double value, const double sigma,
                                 const APLCON::Variable_Settings_t& settings)
{
  if(sigma == 0) {
    throw logic_error("Measured variables need non-zero sigma. By definition, they are unmeasured then.");
  }
  if(settings.StepSize == 0) {
    throw logic_error("Measured variables need non-zero step size. By definition, they are fixed then.");
  }
  AddVariable(name, value, sigma, settings);
}

void APLCON::AddUnmeasuredVariable(const string &name, const double value,
                                   const APLCON::Variable_Settings_t& settings)
{
  if(settings.StepSize == 0) {
    throw logic_error("Unmeasured variables need non-zero step size. By definition, they are fixed then.");
  }
  // unmeasured variables have a sigma of 0
  AddVariable(name, value, 0, settings);
}

void APLCON::AddFixedVariable(const string &name, const double value, const double sigma,
                              const Distribution_t &distribution)
{
  if(sigma == 0) {
    throw logic_error("Fixed variables need non-zero sigma. By definition, they are unmeasured then.");
  }
  // fixed variables have stepSize of 0
  // and limits don't apply (probably?)
  AddVariable(name, value, sigma, {distribution, Variable_Settings_t::Default.Limit, 0});
}

void APLCON::AddVariable(const string &name, const double value, const double sigma,
                         const APLCON::Variable_Settings_t& settings)
{
  // check if variable already exists
  CheckMapKey("Variable", name, variables);


  variable_t var;
  var.StoredValues.emplace_back(value);
  var.StoredSigmas.emplace_back(sigma);
  var.Settings.emplace_back(settings);

  variables[name] = var;
  initialized = false;
}

void APLCON::LinkVariable(const string &name,
                          const std::vector<double*> &values,
                          const std::vector<double*> &sigmas,
                          const std::vector<APLCON::Variable_Settings_t> &settings) {
  // linked sigmas are easy
  auto it = LinkVariable(name, values, sigmas.size(), settings);
  it->second.Sigmas = sigmas;
}

void APLCON::LinkVariable(const string &name,
                          const std::vector<double*>& values,
                          const std::vector<double>& sigmas,
                          const std::vector<APLCON::Variable_Settings_t> &settings) {
  // internally stored sigmas
  // be nice and offer to set sigma for all values to one single provided value
  vector<double> sigmas_(1);
  if(sigmas.size()==1) {
    sigmas_.resize(values.size(), sigmas[0]);
  }
  else {
    sigmas_ = sigmas; // copy
  }
  auto it = LinkVariable(name, values, sigmas_.size(), settings);
  it->second.StoredSigmas = sigmas_;
}

APLCON::variables_t::iterator APLCON::LinkVariable(const string &name,
                          const std::vector<double*>& values,
                          const size_t sigmas_size,
                          const std::vector<APLCON::Variable_Settings_t> &settings) {

  CheckMapKey("Linked Variable", name, variables);

  const size_t n = values.size();
  if(n==0) {
    throw std::logic_error("At least one value should be linked");
  }

  variable_t var;
  if(sigmas_size != n) {
    throw std::logic_error("Sigmas size does not match number of provided values");
  }

  // check settings, set to defaults if none given
  if(settings.size() == 0) {
    var.Settings.resize(n, Variable_Settings_t::Default);
  }
  else if(settings.size() == 1) {
    var.Settings.resize(n, settings[0]);
  }
  else if(settings.size() != n) {
    // 0, 1 or it must fit to values...
    throw std::logic_error("Settings size does not match number of provided values");
  }
  else {
    var.Settings = settings;
  }

  // set values, that's easy
  var.Values = values;

  auto p = variables.insert(make_pair(name, var));
  return p.first;
}

void APLCON::SetCovariance(const string &var1, const string &var2, const double cov)
{
  if(var1.empty() || var2.empty()) {
    throw logic_error("Variable names cannot be empty strings");
  }
  // assume covariance is gonna changed
  initialized = false;

  // search the covariances map
  // but note that the pairs are symmetric (since the covariance matrix is)
  const pair<string, string>& p1 = make_pair(var1, var2);
  const pair<string, string>& p2 = make_pair(var2, var1);
  auto it1 = covariances.find(p1);
  if(it1 != covariances.end()) {
    //it1->second = cov;
    return;
  }
  auto it2 = covariances.find(p2);
  if(it2 != covariances.end()) {
    //it2->second = cov;
    return;
  }
  // not found, then add it
  //covariances.insert(make_pair(p1, cov));
}



APLCON::Result_t APLCON::DoFit()
{
  // ensure that APLCON is properly initialized
  // and the value vectors X, F, V
  Init();

  // the main convergence loop
  int aplcon_ret = -1;
  do {
    // evaluate the constraints F_func and
    // store results in F via iterator F_it
    auto F_it = F.begin();
    for(size_t i=0; i<F_func.size(); ++i) {
      for(const auto& v : F_func[i]()) {
        *F_it = v;
        ++F_it;
      }
    }
    // call APLCON iteration
    c_aplcon_aploop(X.data(), V.data(), F.data(), &aplcon_ret);
  }
  while(aplcon_ret<0);

  // now, after the loop, retrieve the results
  Result_t result = Result_t::Default;

  // make some evil static_cast, but it's way shorter than switch statement
  if(aplcon_ret >= static_cast<int>(Result_Status_t::_Unknown)) {
    throw logic_error("Unkown return value after APLCON fit");
  }
  result.Status = static_cast<Result_Status_t>(aplcon_ret);

  // return default result if not successful
  if(result.Status != Result_Status_t::Success)
    return result;

  // now retrieve everything from APLCON,
  // we don't do any calculation on our own
  // since we are a stupid wrapper :)

  // retrieve some info about the fit (directly copy to struct field if possible)
  float chi2, pval;
  // chndpv and apstat both return the resulting chi2,
  // but the latter returns it with double precision
  c_aplcon_chndpv(&chi2,&result.NDoF,&pval);
  result.Probability = pval;
  c_aplcon_apstat(&result.ChiSquare, &result.NFunctionCalls, &result.NIterations);

  // get the pulls from APLCON
  vector<double> pulls(X.size());
  c_aplcon_appull(pulls.data());

  // now we're ready to fill the result.Variables vector
  // we fill it in the same order as the X vector
  // which makes debug output from  APLCON comparable to dumps of this structure
  result.Variables.reserve(X.size());

  for(const auto& it_var : variables) {
    const string& name = it_var.first;
    const variable_t& before = it_var.second;

    for(size_t k=0;k<before.Values.size();k++) {
      size_t i = before.XOffset+k;

      Result_Variable_t var;
      stringstream s_name;
      s_name << name;
      if(before.Values.size()>1) {
        s_name << "[" << k << "]";
      }
      var.Name = s_name.str();
      var.Value = {*(before.Values[k]), X[i]};

      const size_t V_i = (i+1)*(i+2)/2-1;
      // pow/sqrt of covariance is actually the only calcution the wrapper does
      // the rest is done by APLCON...
      const double after_sigma = sqrt(V[V_i]);
      var.Sigma = {*(before.Sigmas[k]), after_sigma};

      // only copy stuff back if variable is not internally stored
      // which is indicated by an empty store
      if(before.StoredValues.empty())
        *(before.Values[k]) = X[i];
      if(before.StoredSigmas.empty())
        *(before.Sigmas[k]) = after_sigma;

      // copy the covariances, respecting that V is symmetrized
      var.Covariances.Before.resize(X.size());
      var.Covariances.After.resize(X.size());
      for(size_t j=0; j<X.size(); j++) {
        const size_t i_ = i<j ? i : j;
        const size_t j_ = i>j ? i : j;
        const size_t V_i = j_*(j_+1)/2 + i_;
        var.Covariances.Before[j] = V_before[V_i];
        var.Covariances.After[j]  = V[V_i];
      }
      var.Pull = pulls[i];
      var.Settings = before.Settings[k];

      // iterating over variables should be the right order
      result.Variables.push_back(move(var));
    }
  }

  // copy just the names of the constraints
  result.Constraints.reserve(constraints.size());
  for(const auto& c : constraints) {
    Result_Constraint_t con;
    con.Name = c.first;
    con.Number = c.second.Number;
    result.Constraints.emplace_back(con);
  }
  result.Name = instance_name;
  return result;
}


void APLCON::Init()
{
  if(initialized && instance_id == instance_lastfit) {
    // reset APLCON for next fit
    c_aplcon_aplcon(nVariables, nConstraints);

    // copy again the linked variables to X,
    for(auto& it_var : variables) {
      const variable_t& var = it_var.second;
      auto X_offset = X.begin() + var.XOffset;
      auto dereference = [] (double* d) {return *d;};
      transform(var.Values.begin(),var.Values.end(), X_offset, dereference);
    }

    // TODO: and sigmas to diagonal V after original V_before?
    V = V_before;

    return;
  }


  // build the storage arrays for APLCON

  // X are simply the start values, but also track the
  // map of variables names to index in X (as offsets)
  // this is used to create the double pointer arrays
  // for the constraints later and also to unmap results of APLCON in DoFit
  nVariables = 0;
  X.clear();
  V.clear();
  X.reserve(2*variables.size()); // best guess, nVariables will be known later
  V.reserve((X.size() << 1)/2); // another estimate for V's size N^2/2 =~ N*(N+1)/2
  for(auto& it_var : variables) {
    variable_t& var = it_var.second;
    size_t offset = X.size();
    var.XOffset = offset;

    // in Init is the best place to do this pointer business,
    // so make it for the internally stored variables here
    // we need to work on Values/Sigmas independently, since Sigmas are only optionally linked
    APLCON_::copy_pointers(var.StoredValues, var.Values);
    APLCON_::copy_pointers(var.StoredSigmas, var.Sigmas);

    // now, externally linked variables and internally stored can
    // be treated equally
    nVariables += var.Values.size();
    for(size_t i=0;i<var.Values.size();i++) {
      X.push_back(*(var.Values[i])); // copy initial X values

      // take care of diagonal elements in V
      // set off diagonal to zero
      const size_t j = offset+i;
      const size_t V_j = (j+1)*(j+2)/2-1;
      V.resize(V_j+1, 0);
      V[V_j] = pow(*(var.Sigmas[i]),2);
    }
  }



  // F will be set by APLCON iteration loop in DoFit
  // F_func are bound to the double pointers which we know
  // since X is now finally allocated in memory
  nConstraints = 0;
  F_func.clear();
  F_func.reserve(constraints.size());
  for(auto& c_map : constraints) {
    // build the vector of double pointers
    constraint_t& constraint = c_map.second;
    vector< vector<const double*> > args;
    args.reserve(constraint.VariableNames.size()); // args usually smaller, but probably not larger (but not excluded)
    for(const string& varname : constraint.VariableNames) {
      const auto& index = variables.find(varname);
      if(index == variables.end()) {
        throw logic_error("Constraint '"+c_map.first+"' refers to unknown variable '"+varname+"'");
      }
      const variable_t& var = index->second;
      // check if constraint fits to variables
      if(constraint.WantsDouble && var.Values.size()>1) {
        stringstream msg;
        msg << "Constraint '" << c_map.first << "' wants only single double arguments, "
            << "but '" << varname << "' consists of " << var.Values.size() << " (i.e. more than 1) values.";
        throw logic_error(msg.str());
      }
      // build the vector of pointers to X values
      vector<const double*> p(var.Values.size());
      const auto& X_offset = X.begin()+var.XOffset;
      transform(X_offset, X_offset+p.size(), p.begin(), APLCON_::make_pointer<double>);
      // and store it in args
      args.push_back(p);
    }
    const auto& func = bind(constraint.Function, args);
    // now, since we have bound the func, we can execute it once
    // to determine the returned number of values and
    // thus obtain the number of constraints
    const vector<double>&  r = func();
    nConstraints += r.size();
    constraint.Number = r.size();
    F_func.push_back(func);
  }
  F.resize(nConstraints);

  //  // V filled with off-diagonal elements from covariances
  //  for(const auto& c_map : covariances) {
  //    const pair<string, string>& vars = c_map.first;
  //    auto i1 = X_s2i.find(vars.first);
  //    if(i1 == X_s2i.end()) {
  //      throw logic_error("Covariance variable '"+vars.first+"' not found");
  //    }
  //    auto i2 = X_s2i.find(vars.second);
  //    if(i2 == X_s2i.end()) {
  //      throw logic_error("Covariance variable '"+vars.second+"' not found");
  //    }
  //    if(i1==i2) {
  //      throw logic_error("Use variable's Sigma field to define uncertainties.");
  //    }
  //    int i = i1->second;
  //    int j = i2->second;
  //    if(i>j) {
  //      swap(i,j);
  //    }
  //    const size_t V_offset = j*(j+1)/2;
  //    V[V_offset+i] = c_map.second;
  //  }

  // finally we know the number of variables and constraints
  // so we can setup APLCON itself
  c_aplcon_aplcon(nVariables, nConstraints);

  c_aplcon_aprint(0,fit_settings.DebugLevel); // default output on LUNP 0
  if(isfinite(fit_settings.ConstraintAccuracy))
    c_aplcon_apdeps(fit_settings.ConstraintAccuracy);
  if(fit_settings.MaxIterations>=0)
    c_aplcon_apiter(fit_settings.MaxIterations);
  if(isfinite(fit_settings.MeasuredStepSizeFactor))
    c_aplcon_apderf(fit_settings.MeasuredStepSizeFactor);
  if(isfinite(fit_settings.UnmeasuredStepSizeFactor))
    c_aplcon_apderu(fit_settings.UnmeasuredStepSizeFactor);
  if(isfinite(fit_settings.MinimalStepSizeFactor))
    c_aplcon_apdlow(fit_settings.MinimalStepSizeFactor);

  for(const auto& it_var : variables) {
    const variable_t& var = it_var.second;
    for(size_t j=0;j<var.Settings.size();j++) {
      const Variable_Settings_t& s = var.Settings[j];
      int i = j+var.XOffset;
      // setup APLCON variable specific things
      switch (s.Distribution) {
      case APLCON::Distribution_t::Gaussian:
        // thats the APLCON default, nothing must be called
        break;
      case APLCON::Distribution_t::Poissonian:
        c_aplcon_apoiss(i);
        break;
      case APLCON::Distribution_t::LogNormal:
        c_aplcon_aplogn(i);
        break;
      case APLCON::Distribution_t::SquareRoot:
        c_aplcon_apsqrt(i);
        break;
        // APLCON exposes even more transformations (see wrapper),
        // but they're not mentioned in the README...
      default:
        break;
      }

      if(isfinite(s.Limit.Low) && isfinite(s.Limit.High))
        c_aplcon_aplimt(i, s.Limit.Low, s.Limit.High);
      if(isfinite(s.StepSize))
        c_aplcon_apstep(i, s.StepSize);
    }
  }

  // save a copy for later
  V_before = V;

  // remember that this instance has inited APLCON
  initialized = true;
  instance_lastfit = instance_id;
}



// finally the ostream implementation for nice (debug) printout

