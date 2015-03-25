#include "APLCON.hpp"

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <utility>
#include <sstream>
#include <iomanip>


extern "C" {
#include "wrapper/APLCON.h"
}

using namespace std;

int APLCON::instance_counter = 0;
int APLCON::instance_lastfit = 0;
const APLCON::Limit_t APLCON::Limit_t::NoLimit = {
  numeric_limits<double>::quiet_NaN(),
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
                                 const APLCON::Distribution_t distribution,
                                 const Limit_t limit,
                                 const double stepSize)
{
  if(sigma == 0) {
    throw logic_error("Measured variables need non-zero sigma. By definition, they are unmeasured then.");
  }
  if(stepSize == 0) {
    throw logic_error("Measured variables need non-zero step size. By definition, they are fixed then.");
  }
  AddVariable(name, value, sigma, distribution, limit, stepSize);
}

void APLCON::AddUnmeasuredVariable(const string &name, const double value,
                                   const Limit_t limit,
                                   const double stepSize)
{
  if(stepSize == 0) {
    throw logic_error("Unmeasured variables need non-zero step size. By definition, they are fixed then.");
  }
  // unmeasured variables have a sigma of 0
  AddVariable(name, value, 0, Distribution_t::Gaussian, limit, stepSize);
}

void APLCON::AddFixedVariable(const string &name, const double value, const double sigma,
                              const Distribution_t distribution)
{
  if(sigma == 0) {
    throw logic_error("Fixed variables need non-zero sigma. By definition, they are unmeasured then.");
  }
  // fixed variables have stepSize of 0
  // and limits don't apply (probably?)
  AddVariable(name, value, sigma, distribution, Limit_t::NoLimit, 0);
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
    // evaluate the constraints
    for(size_t i=0; i<F.size(); ++i) {
      F[i] = F_func[i]();
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
  result.Variables.reserve(variables.size());

  auto it_vars = variables.cbegin();
  for(size_t i=0; i<X.size(); ++i, ++it_vars) {
    const string& name = it_vars->first;
    const Variable_t& before = it_vars->second;

    Result_Variable_t var;
    var.Name = name;
    var.Value = {before.Value, X[i]};
    const size_t V_i = (i+1)*(i+2)/2-1;
    var.Sigma = {before.Sigma, sqrt(V[V_i])};
    //
    var.Covariances.Before.resize(X.size());
    var.Covariances.After.resize(X.size());
    for(size_t j=0; j<X.size(); j++) {
      const size_t i_ = i<j ? i : j;
      const size_t j_ = i>j ? i : j;
      const size_t V_i = j_*(j_+1)/2 + i_;
      var.Covariances.Before[j] = V_before[V_i];
      var.Covariances.After[j] = V[V_i];
    }
    var.Pull = pulls[i];
    var.Settings = before.Settings;

    // iterating over variables should be the right order
    result.Variables.push_back(var);
  }

  result.Constraints.reserve(constraints.size());
  for(const auto& c : constraints) {
    result.Constraints.emplace_back(c.first);
  }
  result.Name = instance_name;
  return result;
}

void APLCON::Init()
{
  // tell APLCON the number of variables and the number of constraints
  // this call is always needed before the c_aplcon_aploop calls
  c_aplcon_aplcon(variables.size(), constraints.size());

  if(initialized && instance_id == instance_lastfit)
    return;

  // setup APLCON itself
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

  // build the storage arrays for APLCON

  // prepare the symmetric covariance matrix
  const size_t n = variables.size();
  V.resize(n*(n+1)/2);
  fill(V.begin(), V.end(), 0);

  // X are simpy the start values, but also track the
  // map of variables names to index in X and vice-versa
  // this is used to create the double pointer arrays
  // for the constraints later and also to remap results of APLCON in DoFit
  X.resize(variables.size());
  auto it_vars = variables.cbegin(); // use const iterator
  X_s2i.clear(); // clear the map
  for(size_t i = 0; i < X.size() && it_vars != variables.cend() ; ++i, ++it_vars) {
    const Variable_t& var = it_vars->second;
    X[i] = var.Value; // copy initial values
    X_s2i[it_vars->first] = i;
    const size_t V_i = (i+1)*(i+2)/2-1;
    V[V_i] = pow(var.Sigma,2);

    // setup APLCON variable specific things
    switch (var.Settings.Distribution) {
    case Distribution_t::Gaussian:
      // thats the APLCON default
      break;
    case Distribution_t::Poissonian:
      c_aplcon_apoiss(i);
      break;
    case Distribution_t::LogNormal:
      c_aplcon_aplogn(i);
      break;
    case Distribution_t::SquareRoot:
      c_aplcon_apsqrt(i);
      break;
    // APLCON exposes even more transformations (see wrapper),
    // but they're not mentioned in the README...
    default:
      break;
    }

    if(isfinite(var.Settings.Limit.Low) && isfinite(var.Settings.Limit.High))
      c_aplcon_aplimt(i, var.Settings.Limit.Low, var.Settings.Limit.High);
    if(isfinite(var.Settings.StepSize))
      c_aplcon_apstep(i, var.Settings.StepSize);
  }

  // F will be set by APLCON iteration loop in DoFit
  // F_func are bound to the double pointers which we know
  // since X is now finally allocated in memory
  F.resize(constraints.size());
  F_func.clear();
  F_func.reserve(constraints.size());
  for(const auto& c_map : constraints) {
    // build the vector of double pointers
    const constraint_t& c = c_map.second;
    vector<const double*> args;
    args.reserve(c.VariableNames.size());
    for(const string& varname : c.VariableNames) {
      auto index = X_s2i.find(varname);
      if(index == X_s2i.end()) {
        throw logic_error("Constraint '"+c_map.first+"' refers to unknown variable '"+varname+"'");
      }
      args.emplace_back(&X[index->second]);
    }
    F_func.emplace_back(bind(c.Function, args));
  }

  // V filled with off-diagonal elements from covariances
  for(const auto& c_map : covariances) {
    const pair<string, string>& vars = c_map.first;
    auto i1 = X_s2i.find(vars.first);
    if(i1 == X_s2i.end()) {
      throw logic_error("Covariance variable '"+vars.first+"' not found");
    }
    auto i2 = X_s2i.find(vars.second);
    if(i2 == X_s2i.end()) {
      throw logic_error("Covariance variable '"+vars.second+"' not found");
    }
    if(i1==i2) {
      throw logic_error("Use variable's Sigma field to define uncertainties.");
    }
    int i = i1->second;
    int j = i2->second;
    if(i>j) {
      swap(i,j);
    }
    const size_t offset = j*(j+1)/2;
    V[offset+i] = c_map.second;
  }
  // save a copy for later
  V_before = V;

  // remember that this instance has inited APLCON
  initialized = true;
  instance_lastfit = instance_id;
}

void APLCON::AddVariable(const string &name, const double value, const double sigma,
                         const APLCON::Distribution_t distribution,
                         const Limit_t limit,
                         const double stepSize)
{
  // check if variable already exists
  TestName("Variable", name, variables);

  Variable_t var;
  var.Value = value;
  var.Sigma = sigma;
  var.Settings.Distribution = distribution;
  var.Settings.Limit = limit;
  var.Settings.StepSize = stepSize;

  // add the variable to the map
  variables[name] = var;

  initialized = false;
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
