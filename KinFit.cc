#include "KinFit.h"

#include <stdexcept>
#include <cmath>
#include <iostream>

extern "C" {
#include "APLCON_wrapper.h"
}


using namespace std;

int KinFit::instance_counter = 0;
int KinFit::instance_lastfit = 0;


void KinFit::AddMeasuredVariable(const std::string &name, const double value, const double sigma,
                                 const KinFit::Distribution_t distribution,
                                 const double lowerLimit, const double upperLimit,
                                 const double stepSize)
{
  if(sigma == 0) {
    throw logic_error("Measured variables need non-zero sigma. By definition, they are unmeasured then.");
  }
  if(stepSize == 0) {
    throw logic_error("Measured variables need non-zero step size. By definition, they are fixed then.");
  }
  AddVariable(name, value, sigma, distribution, lowerLimit, upperLimit, stepSize);
}

void KinFit::AddUnmeasuredVariable(const string &name, const double value,
                                   const double lowerLimit, const double upperLimit,
                                   const double stepSize)
{
  if(stepSize == 0) {
    throw logic_error("Unmeasured variables need non-zero step size. By definition, they are fixed then.");
  }
  // unmeasured variables have a sigma of 0
  AddVariable(name, value, 0, Distribution_t::Gaussian, lowerLimit, upperLimit, stepSize);
}

void KinFit::AddFixedVariable(const string &name, const double value, const double sigma,
                              const Distribution_t distribution)
{
  if(sigma == 0) {
    throw logic_error("Fixed variables need non-zero sigma. By definition, they are unmeasured then.");
  }
  // fixed variables have stepSize of 0
  // and limits don't apply
  AddVariable(name, value, sigma, distribution,
              std::numeric_limits<double>::quiet_NaN(),
              std::numeric_limits<double>::quiet_NaN(),
              0
              );
}

KinFit::Result_t KinFit::DoFit()
{
  // ensure that APLCON is properly initialized
  // and the value vectors X, F, V
  Init();

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

}

void KinFit::Init()
{
  if(initialized && instance_id == instance_lastfit)
    return;

  // tell APLCON the number of variables and the number of constraints
  c_aplcon_aplcon(variables.size(), constraints.size());

  // TODO: init more APLCON stuff like step sizes and so on...

  // build the storage arrays for APLCON

  // very easy, copy the covariance matrix
  V = covariances;

  // during init for X, also track the
  // map of variables names to index in X
  // this is used to create the double pointer arrays
  // for the constraints later
  X.resize(variables.size());
  auto it_vars = variables.cbegin(); // use const iterator
  map<string, size_t> X_indices;
  for(size_t i = 0; i < X.size() && it_vars != variables.cend() ; ++i, ++it_vars) {
    X[i] = it_vars->second; // copy initial values
    X_indices[it_vars->first] = i;
  }

  // F will be set by APLCON iteration loop in DoFit
  // F_func are bound to the double pointers which we know
  // since X is now finally allocated in memory
  F.resize(constraints.size());
  F_func.reserve(constraints.size());
  for(const auto& c_map : constraints) {
    // build the vector of double pointers
    const constraint_t& c = c_map.second;
    vector<const double*> args;
    args.reserve(c.VariableNames.size());
    for(const string& varname : c.VariableNames) {
      auto index = X_indices.find(varname);
      if(index == X_indices.end()) {
        throw logic_error("Constraint '"+c_map.first+"' refers to unknown variable '"+varname+"'");
      }
      args.emplace_back(&X[index->second]);
    }
    F_func.emplace_back(bind(c.Function, args));
  }

  // remember that this instance has inited APLCON
  initialized = true;
  instance_lastfit = instance_id;
}

void KinFit::AddVariable(const string &name, const double value, const double sigma,
                         const KinFit::Distribution_t distribution,
                         const double lowerLimit, const double upperLimit,
                         const double stepSize)
{
  // check if variable already exists
  TestName("Variable", name, variables);

  // add the variable to the map
  variables[name] = value;
  const size_t n = variables.size();

  // resize the symmetric covariance matrix and set the sigma
  const size_t cov_size = n*(n+1)/2;
  covariances.resize(cov_size, 0);
  covariances[cov_size-1] = pow(sigma,2);

  // save the other information
  distributions.push_back(distribution);
  limits.push_back(make_pair(lowerLimit, upperLimit));
  stepSizes.push_back(stepSize);

  initialized = false;
}

