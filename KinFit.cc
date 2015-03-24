#include "KinFit.h"

#include <stdexcept>
#include <cmath>

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

/*void KinFit::AddConstraint(const string &name, const constraint_t &constraint)
{
  // check if variable already exists
  if(constraints.find(name) != constraints.end()) {
    throw logic_error("Constraint already added");
  }
  constraints[name] = constraint;
  initialized = false;
}*/

KinFit::Result_t KinFit::DoFit()
{
  // ensure that APLCON is properly initialized
  // and the value vectors X, F, V
  Init();

  int aplcon_ret = -1;
  do {
    // evaluate the constraints
    //vector<double> evals(constraints.size());
    //auto it3 = constraints.begin();
    //auto it4 = evals.begin();
    //c_aplcon_aploop(values.data(), covariances.data(), F.data(), &IRET);
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

