#include "KinFit.h"

#include <stdexcept>
#include <cmath>

extern "C" {
#include "APLCON_wrapper.h"
}


using namespace std;


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
              std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
              0
              );
}

void KinFit::AddVariable(const string &name, const double value, const double sigma,
                         const KinFit::Distribution_t distribution,
                         const double lowerLimit, const double upperLimit,
                         const double stepSize)
{
  // check if variable already exists
  if(variables.find(name) != variables.end()) {
    throw logic_error("Variable already exists");
  }

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
}

