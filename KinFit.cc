#include "KinFit.h"

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <utility>
#include <sstream>
#include <iomanip>


extern "C" {
#include "APLCON_wrapper.h"
}


using namespace std;

int KinFit::instance_counter = 0;
int KinFit::instance_lastfit = 0;
const KinFit::Limit_t KinFit::Limit_t::NoLimit = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};

void KinFit::AddMeasuredVariable(const std::string &name, const double value, const double sigma,
                                 const KinFit::Distribution_t distribution,
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

void KinFit::AddUnmeasuredVariable(const string &name, const double value,
                                   const Limit_t limit,
                                   const double stepSize)
{
  if(stepSize == 0) {
    throw logic_error("Unmeasured variables need non-zero step size. By definition, they are fixed then.");
  }
  // unmeasured variables have a sigma of 0
  AddVariable(name, value, 0, Distribution_t::Gaussian, limit, stepSize);
}

void KinFit::AddFixedVariable(const string &name, const double value, const double sigma,
                              const Distribution_t distribution)
{
  if(sigma == 0) {
    throw logic_error("Fixed variables need non-zero sigma. By definition, they are unmeasured then.");
  }
  // fixed variables have stepSize of 0
  // and limits don't apply
  AddVariable(name, value, sigma, distribution, Limit_t::NoLimit, 0);
}

void KinFit::AddCovariance(const string &var1, const string &var2, const double cov)
{
  if(var1.empty() || var2.empty()) {
    throw logic_error("Variable names cannot be empty strings");
  }
  const pair<string, string>& p1 = make_pair(var1, var2);
  const pair<string, string>& p2 = make_pair(var2, var1);
  if(covariances.find(p1) != covariances.end() ||
     covariances.find(p2) != covariances.end()) {
    throw logic_error("Covariance between '"+var1+"' and '"+var2+"' already added.");
  }
  covariances[p1] = cov;
  initialized = false;
}

KinFit::Result_t KinFit::DoFit()
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

void KinFit::Init()
{
  // tell APLCON the number of variables and the number of constraints
  // this call is always needed before the c_aplcon_aploop calls
  c_aplcon_aplcon(variables.size(), constraints.size());

  if(initialized && instance_id == instance_lastfit)
    return;

  c_aplcon_aprint(0,0); // no output on stdout from now on

  // TODO: init more APLCON stuff like step sizes and so on...

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
  //X_i2s.resize(variables.size());
  for(size_t i = 0; i < X.size() && it_vars != variables.cend() ; ++i, ++it_vars) {
    const Variable_t& var = it_vars->second;
    X[i] = var.Value; // copy initial values
    X_s2i[it_vars->first] = i;
    const size_t V_i = (i+1)*(i+2)/2-1;
    V[V_i] = pow(var.Sigma,2);
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

void KinFit::AddVariable(const string &name, const double value, const double sigma,
                         const KinFit::Distribution_t distribution,
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

std::ostream& operator<< (std::ostream& o, const KinFit::Result_Status_t& s) {
  switch (s) {
  case KinFit::Result_Status_t::Success:
    o << "Success";
    break;
  case KinFit::Result_Status_t::NoConvergence:
    o << "NoConvergence";
    break;
  case KinFit::Result_Status_t::TooManyIterations:
    o << "TooManyIterations";
    break;
  case KinFit::Result_Status_t::UnphysicalValues:
    o << "UnphysicalValues";
    break;
  case KinFit::Result_Status_t::NegativeDoF:
    o << "NegativeDoF";
    break;
  case KinFit::Result_Status_t::OutOfMemory:
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
string stringify_contraints(const vector<KinFit::Result_Variable_t>& variables, F f) {
  const string& in = "   "; // indent
  const int w = 15;
  stringstream o;
  o << in<<in << "Covariances: " << endl;
  o << in<<in << setw(w) << " ";
  for(size_t i=0;i<variables.size();i++) {
    stringstream i_;
    i_ << "(" << i << ")";
    o << setw(w) << i_.str();
  }
  o << endl;
  for(size_t i=0;i<variables.size();i++) {
    const KinFit::Result_Variable_t& v = variables[i];
    stringstream i_;
    i_ << "(" << i << ") ";
    o << in<<in << setw(4) << i_.str() << " " << setw(10) << v.Name;
    const vector<double>& cov = f(v);
    for(size_t j=0;j<cov.size();j++) {
      o << setw(w) << cov[j];
    }
    o << endl;
  }
  o << in<<in << setw(w) << " ";
  for(size_t i=0;i<variables.size();i++) {
    stringstream i_;
    i_ << "(" << i << ")";
    o << setw(w) << i_.str();
  }
  o << endl;
  o << endl;
  return o.str();
}

ostream& operator<< (ostream& o, const KinFit::Result_t& r) {
  const string& in = "   "; // indent
  const string& ma = ">> "; // marker
  const int w = 10;

  // general info
  o << ma << (r.Name==""?"KinFit":r.Name) << " with " << r.Variables.size() << " variables and "
    << r.Constraints.size() << " constraints:" << endl;
  o << in << r.Status << " after " << r.NIterations << " iterations, " << r.NFunctionCalls << " function calls " << endl;
  o << in << "Constraints: " <<
       stringify(r.Constraints, [](const string& v) {return v;}) << endl << endl;

  // print stuff before the Fit
  o << in<<ma << "Before Fit:" << endl;
  o << in<<in << setw(w) << "Name"
    << setw(w) << "Value"
    << setw(w) << "Sigma"
    << endl;
  for(const auto& v : r.Variables) {
    o << in<<in << setw(w) << v.Name
      << setw(w) << v.Value.Before << "  "
      << setw(w) << v.Sigma.Before << "  "
      << endl;
  }
  o << endl;

  o << stringify_contraints(r.Variables, [](const KinFit::Result_Variable_t& v) {return v.Covariances.Before;});

  // print stuff after the fit
  o << in<<ma << "After Fit:" << endl;
  o << in<<in << setw(w) << "Name"
    << setw(w) << "Value"
    << setw(w) << "Sigma"
    << setw(w) << "Pull"
    << endl;
  for(const auto& v : r.Variables) {
    o << in<<in << setw(w) << v.Name
      << setw(w) << v.Value.After << "  "
      << setw(w) << v.Sigma.After << "  "
      << setw(w) << v.Pull
      << endl;
  }
  o << endl;

  o << stringify_contraints(r.Variables, [](const KinFit::Result_Variable_t& v) {return v.Covariances.After;});

  return o;
}

