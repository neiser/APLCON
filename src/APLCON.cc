#include "APLCON.hpp"

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <utility>
#include <sstream>
#include <iomanip>
#include <limits>


extern "C" {
#include "wrapper/APLCON.h"
}

using namespace std;

int APLCON::instance_counter = 0;
int APLCON::instance_lastfit = 0;

const double APLCON::NaN = numeric_limits<double>::quiet_NaN();
const APLCON::Variable_Settings_t APLCON::Variable_Settings_t::Default = {
  APLCON::Distribution_t::Gaussian,
  {
    -numeric_limits<double>::infinity(),
    numeric_limits<double>::infinity()
  },
  numeric_limits<double>::quiet_NaN()
};
// set in Init() method
const APLCON::Fit_Settings_t APLCON::Fit_Settings_t::Default = {
  0,  // no debug printout
  -1, // default max iterations
  numeric_limits<double>::quiet_NaN(),
  numeric_limits<double>::quiet_NaN(),
  numeric_limits<double>::quiet_NaN(),
  numeric_limits<double>::quiet_NaN()
};

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


  Variable_t var;
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
  CheckMapKey("Linked Variable", name, variables);

  const size_t n = values.size();
  if(n==0) {
    throw std::logic_error("At least one value should be linked");
  }

  Variable_t var;
  if(sigmas.size() != n) {
    throw std::logic_error("Sigmas size does not match number of of provided values");
  }

  if(settings.size() == 0) {
    var.Settings.resize(n, Variable_Settings_t::Default);
  }
  else if(settings.size() != n) {
    throw std::logic_error("Sigmas size does not match number of provided values");
  }
  else {
    var.Settings = settings;
  }

  // set values and sigmas
  var.Values = values;
  var.Sigmas = sigmas;

  variables.insert(make_pair(name, var));
}

void APLCON::SetCovariance(const string &var1, const string &var2, const double cov)
{
  if(var1.empty() || var2.empty()) {
    throw logic_error("Variable names cannot be empty strings");
  }
  // assume covariance is gonna changed...
  initialized = false;

  // search the covariances map
  // but note that the pairs are symmetric
  const pair<string, string>& p1 = make_pair(var1, var2);
  const pair<string, string>& p2 = make_pair(var2, var1);
  auto it1 = covariances.find(p1);
  if(it1 != covariances.end()) {
    it1->second = cov;
    return;
  }
  auto it2 = covariances.find(p2);
  if(it2 != covariances.end()) {
    it2->second = cov;
    return;
  }
  // not found, then add it
  covariances.insert(make_pair(p1, cov));
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
  Result_t result;

  // interpret status number according to README
  bool converged = false;
  switch (aplcon_ret) {
  case 0:
    result.Status = Result_Status_t::Success;
    converged = true;
    break;
  case 1:
    result.Status = Result_Status_t::NoConvergence;
    break;
  case 2:
    result.Status = Result_Status_t::TooManyIterations;
    break;
  case 3:
    result.Status = Result_Status_t::UnphysicalValues;
    break;
  case 4:
    result.Status = Result_Status_t::NegativeDoF;
    break;
  case 5:
    result.Status = Result_Status_t::OutOfMemory;
    break;
  default:
    throw logic_error("Unkown return value after APLCON fit");
  }

  // return default result if not successful
  if(!converged)
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

  for(const auto& var_pair : variables) {
    const string& name = var_pair.first;
    const Variable_t& before = var_pair.second;

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
      *(before.Values[k]) = X[i];

      const size_t V_i = (i+1)*(i+2)/2-1;
      var.Sigma = {*(before.Sigmas[k]), sqrt(V[V_i])};

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
    result.Constraints.emplace_back(c.first);
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
    // and sigmas to diagonal V after original V_before
    for(auto& it_var : variables) {
      const Variable_t& var = it_var.second;
      auto X_offset = X.begin() + var.XOffset;
      auto dereference = [] (double* d) {return *d;};
      transform(var.Values.begin(),var.Values.end(), X_offset, dereference);
      V = V_before;
      //transform(var.Sigmas.begin(),var.Sigmas.end(), X_offset, dereference);
    }
    return;
  }


  // build the storage arrays for APLCON

  // X are simpy the start values, but also track the
  // map of variables names to index in X (as offsets)
  // this is used to create the double pointer arrays
  // for the constraints later and also to unmap results of APLCON in DoFit
  nVariables = 0;
  X.clear();
  V.clear();
  X.reserve(2*variables.size()); // best guess, nVariables will be known later
  V.reserve((X.size() << 1)/2); // another estimate for V's size N^2/2 =~ N*(N+1)/2
  for(auto& it_var : variables) {
    Variable_t& var = it_var.second;
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
  for(const auto& c_map : constraints) {
    // build the vector of double pointers
    const constraint_t& c = c_map.second;
    vector< vector<const double*> > args;
    args.reserve(c.VariableNames.size()); // args might be smaller, but probably not larger
    for(const string& varname : c.VariableNames) {
      const auto& index = variables.find(varname);
      if(index == variables.end()) {
        throw logic_error("Constraint '"+c_map.first+"' refers to unknown variable '"+varname+"'");
      }
      const Variable_t& var = index->second;
      // build the vector of pointers to X values
      vector<const double*> p(var.Values.size());
      const auto& X_offset = X.begin()+var.XOffset;
      transform(X_offset, X_offset+p.size(), p.begin(), APLCON_::make_pointer<double>);
      // and store it in args
      args.push_back(p);
    }
    const auto& func = bind(c.Function, args);
    // now, since we have bound the func, we can execute it once
    // to determine the returned number of values and
    // thus obtain the number of constraints
    const vector<double>&  r = func();
    nConstraints += r.size();
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
    const Variable_t& var = it_var.second;
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

const string APLCON::PrintFormatting::Indent = "   ";
const string APLCON::PrintFormatting::Marker = ">> ";
const int APLCON::PrintFormatting::Width = 13;

std::ostream& operator<< (std::ostream& o, const APLCON::Limit_t& l) {
  return o << "(" << l.Low << ", " << l.High << ")";
}

std::ostream& operator<< (std::ostream& o, const APLCON::Distribution_t& d) {
  switch (d) {
  case APLCON::Distribution_t::Gaussian:
    o << "Gaussian";
    break;
  case APLCON::Distribution_t::LogNormal:
    o << "LogNormal";
    break;
  case APLCON::Distribution_t::Poissonian:
    o << "Poissonian";
    break;
  case APLCON::Distribution_t::SquareRoot:
    o << "SquareRoot";
    break;
  default:
    throw logic_error("Unkown Distribution_t in ostream");
    break;
  }
  return o;
}

std::ostream& operator<< (std::ostream& o, const APLCON::Result_Status_t& s) {
  switch (s) {
  case APLCON::Result_Status_t::Success:
    o << "Success";
    break;
  case APLCON::Result_Status_t::NoConvergence:
    o << "NoConvergence";
    break;
  case APLCON::Result_Status_t::TooManyIterations:
    o << "TooManyIterations";
    break;
  case APLCON::Result_Status_t::UnphysicalValues:
    o << "UnphysicalValues";
    break;
  case APLCON::Result_Status_t::NegativeDoF:
    o << "NegativeDoF";
    break;
  case APLCON::Result_Status_t::OutOfMemory:
    o << "OutOfMemory";
    break;
  default:
    throw logic_error("Unkown Result_Status_t in ostream");
    break;
  }

  return o;
}

template<typename T, typename F>
string stringify(const vector<T>& c, F f) {
  stringstream o;
  for(size_t i=0 ; i<c.size() ; i++) {
    o << f(c[i]);
    if(i != c.size()-1)
      o << ", ";
  }
  return o.str();
}

template<typename F>
string stringify_contraints(const vector<APLCON::Result_Variable_t>& variables, const string& in, F f) {
  const int w = APLCON::PrintFormatting::Width;
  const int w_varname = APLCON::PrintFormatting::Width+5;
  stringstream o;
  o << in << "Covariances: " << endl;
  o << in << setw(w_varname) << " ";
  for(size_t i=0;i<variables.size();i++) {
    stringstream i_;
    i_ << "(" << i << ")";
    o << setw(w) << i_.str();
  }
  o << endl;
  for(size_t i=0;i<variables.size();i++) {
    const APLCON::Result_Variable_t& v = variables[i];
    stringstream i_;
    i_ << "(" << i << ") ";
    o << in << setw(4) << i_.str() << " " << setw(w_varname-5) << v.Name;
    const vector<double>& cov = f(v);
    for(size_t j=0;j<cov.size();j++) {
      o << setw(w) << cov[j];
    }
    o << endl;
  }
  o << in << setw(w_varname) << " ";
  for(size_t i=0;i<variables.size();i++) {
    stringstream i_;
    i_ << "(" << i << ")";
    o << setw(w) << i_.str();
  }
  o << endl;
  return o.str();
}

string stringify_variables(const vector<APLCON::Result_Variable_t>& variables, const string& extra_indent = "") {
  stringstream o;
  const int w = APLCON::PrintFormatting::Width;
  const string& in = extra_indent + APLCON::PrintFormatting::Indent;
  const string& ma = extra_indent + APLCON::PrintFormatting::Marker;
  // print stuff before the Fit
  o << ma << "Before Fit:" << endl << endl;
  o << in
    << left << setw(w) << "Name" << right
    << setw(w)   << "Value"
    << setw(w)   << "Sigma"
    << left << setw(2*w) << "   Settings" << right
    << endl;
  for(const APLCON::Result_Variable_t& v : variables) {
    stringstream settings;
    settings << "   " << v.Settings.Distribution << " " << v.Settings.Limit << " " << v.Settings.StepSize;
    o << in
      << left << setw(w) << v.Name << right
      << setw(w) << v.Value.Before
      << setw(w) << v.Sigma.Before
      << left << setw(2*w) << settings.str() << right
      << endl;
  }
  o << endl;

  o << stringify_contraints(variables, in, [](const APLCON::Result_Variable_t& v) {return v.Covariances.Before;});

  // print stuff after the fit
  o << ma << "After Fit:" << endl << endl;
  o << in
    << left << setw(w) << "Variable" << right
    << setw(w) << "Value"
    << setw(w) << "Sigma"
    << setw(w) << "Pull"
    << endl;
  for(const APLCON::Result_Variable_t& v : variables) {
    o << in
      << left << setw(w) << v.Name << right
      << setw(w) << v.Value.After
      << setw(w) << v.Sigma.After
      << setw(w) << v.Pull
      << endl;
  }
  o << endl;

  o << stringify_contraints(variables, in, [](const APLCON::Result_Variable_t& v) {return v.Covariances.After;});

  return o.str();

}

ostream& operator<< (ostream& o, const APLCON::Result_t& r) {
  const string& in = APLCON::PrintFormatting::Indent;
  const string& ma = APLCON::PrintFormatting::Marker;

  // general info
  o << ma << (r.Name==""?"APLCON":r.Name) << " with " << r.Variables.size() << " variables and "
    << r.Constraints.size() << " constraints:" << endl;
  o << in << r.Status << " after " << r.NIterations << " iterations, " << r.NFunctionCalls << " function calls " << endl;
  o << in << "Chi^2 / DoF = " << r.ChiSquare << " / " << r.NDoF << " = " << r.ChiSquare/r.NDoF << endl;
  o << in << "Probability = " << r.Probability << endl;
  o << in << "Constraints: " <<
       stringify(r.Constraints, [](const string& v) {return v;}) << endl << endl;

  if(r.Status == APLCON::Result_Status_t::Success)
    o << stringify_variables(r.Variables, in);
  return o;
}
