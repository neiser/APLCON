#include "KinFit.h"

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <utility>

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
    var.Pull = pulls[i];
    var.Settings = before.Settings;

    // iterating over variables should be the right order
    result.Variables.push_back(var);
  }

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
    size_t offset = j*(j+1)/2;
    V[offset+i] = c_map.second;
  }
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

