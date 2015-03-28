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

// simple AddVariable methods

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

// LinkVariable methods

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
  const size_t sigmas_size = sigmas.size()==1 ? values.size() : sigmas.size();
  auto it = LinkVariable(name, values, sigmas_size, settings);
  it->second.StoredSigmas = sigmas;
  if(sigmas.size()==1) {
    it->second.StoredSigmas.resize(sigmas_size, sigmas[0]);
  }
}

APLCON::variables_t::iterator APLCON::LinkVariable(
    const string &name,
    const std::vector<double*>& values,
    const size_t sigmas_size,
    const std::vector<APLCON::Variable_Settings_t> &settings
    ) {

  CheckMapKey("Linked Variable", name, variables);

  const size_t n = values.size();
  if(n==0) {
    throw std::logic_error("At least one value should be linked");
  }

  if(sigmas_size != n) {
    throw std::logic_error("Sigmas size does not match number of provided values");
  }

  variable_t var;

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

// SetCovariance methods

// assume var1 and var2 are both scalar
void APLCON::SetCovariance(const std::string& var1,
                           const std::string& var2,
                           const double covariance) {

  if(var1 == var2) {
    throw logic_error("Covariance variable names must be different");
  }
  auto it = MakeCovarianceEntry(var1, var2);
  auto& values = it->second.StoredValues;
  values.resize(1);
  values.back() = covariance;
}
void APLCON::SetCovariance(const std::string& var1,
                           const std::string& var2,
                           const std::vector<double>& covariance) {
  auto it = MakeCovarianceEntry(var1, var2);
  if(covariance.empty()) {
    throw logic_error("Empty covariance values given");
  }
  it->second.StoredValues = covariance;
}

void APLCON::SetCovariance(const std::string& var1,
                           const std::string& var2,
                           const std::vector<double*>& covariance) {
  auto it = MakeCovarianceEntry(var1, var2);
  if(covariance.empty()) {
    throw logic_error("Empty covariance pointers given");
  }
  it->second.Values = covariance;
}

APLCON::covariances_t::iterator APLCON::MakeCovarianceEntry(
    const string &var1,
    const string &var2)
{
  if(var1.empty() || var2.empty()) {
    throw logic_error("Covariance variable names cannot be empty strings");
  }
  // we cannot check in general if the variables should have different names
  // since for vector variables, there are still correlations possible
  // finally, this will be checked by Init

  // assume covariance is gonna changed
  // so force Init to execute again
  initialized = false;

  // search the covariances map
  // but note that the pairs are symmetric (since the covariance matrix is)
  // so we search for both possibilities
  // if we find it, we update it...
  const pair<string, string>& p1 = make_pair(var1, var2);
  auto it1 = covariances.find(p1);
  if(it1 != covariances.end()) {
    return it1;
  }
  const pair<string, string>& p2 = make_pair(var2, var1);
  auto it2 = covariances.find(p2);
  if(it2 != covariances.end()) {
    return it2;
  }

  // not found, then add default struct and return
  auto p = covariances.insert(make_pair(p1, covariance_t()));
  return p.first;
}

// Main Fit Routines

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

  for(const auto& it_map : variables) {
    const string& name = it_map.first;
    const variable_t& before = it_map.second;

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

      // pow/sqrt of covariance is actually the only calcution the wrapper does
      // the rest is done by APLCON...
      const size_t V_i = APLCON_::V_ij(i,i);
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
        const size_t V_ij = APLCON_::V_ij(i,j);
        var.Covariances.Before[j] = V_before[V_ij];
        var.Covariances.After[j]  = V[V_ij];
      }

      // calculate the correlations (after building the full covariance matrix),
      // because APLCON does not provide them (only print methods available)
      var.Correlations.Before.resize(X.size());
      var.Correlations.After.resize(X.size());
      for(size_t j=0; j<X.size(); j++) {
        const size_t V_ij = APLCON_::V_ij(i,j);
        const size_t V_ii = APLCON_::V_ij(i,i);
        const size_t V_jj = APLCON_::V_ij(j,j);
        const double prod_before = V_before[V_ii] * V_before[V_jj];
        const double prod = V[V_ii] * V[V_jj];
        var.Correlations.Before[j] = V_before[V_ij]/sqrt(prod_before);
        var.Correlations.After[j]  = V[V_ij]/sqrt(prod);
      }

      // anything else
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
  // check if we can do some quick init
  if(initialized && instance_id == instance_lastfit) {
    // reset APLCON for next fit
    c_aplcon_aplcon(nVariables, nConstraints);

    // reset the covariance matrix
    V = V_before;

    // copy again the linked variables to X
    //
    for(const auto& it_map : variables) {
      const variable_t& var = it_map.second;
      auto X_offset = X.begin() + var.XOffset;
      auto dereference = [] (const double* d) {return *d;};
      transform(var.Values.begin(), var.Values.end(), X_offset, dereference);
      // copy the sigmas to diagonal of V
      APLCON_::transform_to_V(V, var.Sigmas, var.V_ij,
                              [] (double d) {return pow(d,2);});
    }

    // copy again the true non-diagonal covariances
    // the distinction between sigmas and covariances makes the interface hopefully more usable,
    // because always specifying covariances is tedious, but sigmas are crucial
    // for measured variables or, say, constrained fitting
    for(const auto& it_map : covariances) {
      const auto& cov = it_map.second;
      APLCON_::transform_to_V(V, cov.Values, cov.V_ij);
    }
    return;
  }


  // build the storage arrays X, V, F for APLCON

  // X are simply the start values, but also track the
  // map of variables names to index in X (as offsets)
  // this is used to create the double pointer arrays
  // for the constraints later and also to unmap results of APLCON in DoFit
  nVariables = 0;
  X.clear();
  V.clear();
  X.reserve(2*variables.size()); // best guess, nVariables will be known later
  V.reserve((X.size() << 1)/2); // another estimate for V's size N^2/2 =~ N*(N+1)/2
  for(auto& it_map : variables) {
    variable_t& var = it_map.second;
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
    var.V_ij.resize(var.Values.size());
    for(size_t i=0;i<var.Values.size();i++) {
      X.push_back(*(var.Values[i])); // copy initial X values

      // take care of diagonal elements in V,
      // off diagonal elements are set to zero by resize (see also covariance handling later)
      // this is like a push_back but with some padding in between
      const size_t j = offset+i;
      const size_t V_ij = APLCON_::V_ij(j,j);
      var.V_ij[i] = V_ij; // store for later
      V.resize(V_ij+1, 0);
      V.back() = pow(*(var.Sigmas[i]),2);
    }
  }



  // F will be set by APLCON iteration loop in DoFit
  // F_func are bound to the double pointers which we know
  // since X is now finally allocated in memory
  nConstraints = 0;
  F_func.clear();
  F_func.reserve(constraints.size());
  for(auto& it_map : constraints) {
    // build the vector of double pointers
    constraint_t& constraint = it_map.second;
    vector< vector<const double*> > args;
    args.reserve(constraint.VariableNames.size()); // args usually smaller, but probably not larger (but not excluded)
    for(const string& varname : constraint.VariableNames) {
      const auto& it = variables.find(varname);
      if(it == variables.end()) {
        throw logic_error("Constraint '"+it_map.first+"' refers to unknown variable '"+varname+"'");
      }
      const variable_t& var = it->second;
      // check if constraint fits to variables
      if(constraint.WantsDouble && var.Values.size()>1) {
        stringstream msg;
        msg << "Constraint '" << it_map.first << "' wants only single double arguments, "
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

  // V filled with off-diagonal elements from covariances
  // the variables already have their Values pointer correctly filled,
  // so it's size can be used for the variables's dimensionality
  // V is composed of submatrices due to this different dimensionality.
  for(auto& it_map : covariances) {
    const pair<string, string>& varnames = it_map.first;
    covariance_t& cov = it_map.second;

    // create the pointers for internally stored covariances
    APLCON_::copy_pointers(cov.StoredValues, cov.Values);

    const string& cov_name = "<'"+varnames.first+"','"+varnames.second+"'>";

    // the setup of the index mapping filed V_ij in constraint_t
    // differs strongly for off-diagonal vs. diagonal covariance elements

    if(varnames.first == varnames.second) {
      // varnames are equal, only one search required
      auto it = variables.find(varnames.first);
      if(it == variables.end()) {
        throw logic_error("Variable name '"+varnames.first+"' for covariance "+cov_name+" not defined");
      }
      const variable_t& var = it->second;
      const size_t n = var.Values.size();
      if(var.Values.size()==1) {
        throw logic_error("Use sigma to define uncertainty of scalar covariance "+cov_name);
      }
      const size_t v_n = n*(n-1)/2; // expected size of the submatrix without diagonal elements

      if(v_n != cov.Values.size()) {
        stringstream msg;
        msg << "Covariance " << cov_name << " provides " << cov.Values.size()
            << " element" << (cov.Values.size()==1?"":"s") << ", but " << v_n << " covariances needed with"
            << " variable dimensions <" << n << "," << n << ">";
        throw logic_error(msg.str());
      }

      cov.V_ij.reserve(cov.Values.size());
      for(size_t i=0;i<n;i++) {
        for(size_t j=0;j<i;j++) {
          cov.V_ij.push_back(APLCON_::V_ij(i+var.XOffset,j+var.XOffset));
        }
      }
    }
    else {
      // varnames are different, search them both
      auto it1 = variables.find(varnames.first);
      if(it1 == variables.end()) {
        throw logic_error("First variable name '"+varnames.first+"' for covariance "+cov_name+" not defined");
      }

      auto it2 = variables.find(varnames.second);
      if(it2 == variables.end()) {
        throw logic_error("Second variable name '"+varnames.second+"' for covariance "+cov_name+" not defined");
      }
      const variable_t& var1 = it1->second;
      const variable_t& var2 = it2->second;
      const size_t n1 = var1.Values.size();
      const size_t n2 = var2.Values.size();
      const size_t v_n = n1*n2; // expected size of the submatrix, there are no di

      if(v_n != cov.Values.size()) {
        stringstream msg;
        msg << "Covariance " << cov_name << " provides " << cov.Values.size()
            << " element" << (cov.Values.size()==1?"":"s") << ", but " << v_n << " covariances needed with"
            << " variable dimensions <" << n1 << "," << n2 << ">";
        throw logic_error(msg.str());
      }
      // var1 refers to rows, var2 to columns (standard mathematics convention)
      cov.V_ij.reserve(cov.Values.size());
      for(size_t i=0;i<n1;i++) {
        for(size_t j=0;j<n2;j++) {
          cov.V_ij.push_back(APLCON_::V_ij(i+var1.XOffset,j+var2.XOffset));
        }
      }
    }

    // now, with some properly initialized cov.V_ij for each case,
    // we can fill the non-diagonal values of V
    APLCON_::transform_to_V(V, cov.Values, cov.V_ij);

  }

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

  // save a pristine copy for later
  V_before = V;

  // remember that this instance has inited APLCON
  initialized = true;
  instance_lastfit = instance_id;
}
