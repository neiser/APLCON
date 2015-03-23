#ifndef KINFIT_H
#define KINFIT_H

#include <vector>
#include <map>
#include <functional>
#include <limits>

/**
 * @brief The KinFit class
 * Provides a C++11'ish wrapper around
 * V.Blobel's FORTRAN APLCON constrained least squares fitter
 * see http://www.desy.de/~blobel/wwwcondl.html for details
 */
class KinFit
{
public:
  KinFit() {}
  ~KinFit() {}

  enum class Distribution_t {
    Gaussian,
    Poissonian,
    LogNormal,
    SquareRoot
  };

  /**
   * @brief AddMeasuredVariable
   * @param name unique label for variable
   * @param value initial value for variable
   * @param sigma sqrt of diagonal entry in covariance matrix
   * @param distribution optional type of distribution
   * @param lowerLimit lower limit of the variable's value
   * @param upperLimit upper limit of the variable's value
   * @param stepSize step size for numerical derivation
   */
  void AddMeasuredVariable(const std::string& name,
                           const double value = std::numeric_limits<double>::quiet_NaN(),
                           const double sigma = std::numeric_limits<double>::quiet_NaN(),
                           const Distribution_t distribution = Distribution_t::Gaussian,
                           const double lowerLimit = std::numeric_limits<double>::quiet_NaN(),
                           const double upperLimit = std::numeric_limits<double>::quiet_NaN(),
                           const double stepSize = std::numeric_limits<double>::quiet_NaN()
      );
  /**
   * @brief AddUnmeasuredVariable
   * @param name unique label for variable
   * @param value initial value for variable
   * @param lowerLimit lower limit of the variable's value
   * @param upperLimit upper limit of the variable's value
   * @param stepSize step size for numerical derivation
   */
  void AddUnmeasuredVariable(const std::string& name,
                             const double value = std::numeric_limits<double>::quiet_NaN(),
                             const double lowerLimit = std::numeric_limits<double>::quiet_NaN(),
                             const double upperLimit = std::numeric_limits<double>::quiet_NaN(),
                             const double stepSize = std::numeric_limits<double>::quiet_NaN()
      );
  /**
   * @brief AddFixedVariable
   * @param name unique label for variable
   * @param value initial value for variable
   * @param sigma sqrt of diagonal entry in covariance matrix
   * @param distribution optional type of distribution
   */
  void AddFixedVariable(const std::string& name,
                        const double value = std::numeric_limits<double>::quiet_NaN(),
                        const double sigma = std::numeric_limits<double>::quiet_NaN(),
                        const Distribution_t distribution = Distribution_t::Gaussian
      );

  // constraints are functions of the input variables. Should return 0.0 if fulfilled.
  typedef std::function<double (const std::map<std::string, double>&)> constraint_t;

  void AddConstraint(const std::string& name, const constraint_t& constraint);

  void UpdateValues(const std::map<std::string, double>& values);
  void UpdateSigmas(const std::map<std::string, double>& sigmas);

  std::map<std::string, double> GetValues() const;
  std::map<std::string, double> GetSigmas() const;

  enum class Status_t {
    Success,
    NotConverged,
    TooManyIterations,
    UnphysicalValues,
    NegativeDoF,
    OutOfMemory
  };

  struct Result_t {
    Status_t Status;
    double ChiSquare;
    int NDoF;
    double Probability;
    int NIterations;
    int NFunctionCalls;
  };

  Result_t DoFit();

private:
  // values of map represent X (works since map is ordered)
  std::map<std::string, double> variables;
  // track the type of distributions
  std::vector<Distribution_t> distributions;
  // track the limits (NaN if unset)
  std::vector< std::pair<double, double> > limits;
  // represents the symmetric covariance matrix
  std::vector<double> covariances;
   // the constraints as C++11 lambdas
  std::map<std::string, constraint_t> constraints;
  // step sizes for numerical evaluation (zero if fixed, NaN if APLCON default)
  std::vector<double> stepSizes;
  bool initialized;

  void AddVariable(const std::string& name, const double value, const double sigma,
                   const Distribution_t distribution,
                   const double lowerLimit, const double upperLimit,
                   const double stepSize
      );
};

#endif // KINFIT_H
