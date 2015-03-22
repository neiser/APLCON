#ifndef KINFIT_H
#define KINFIT_H

#include <vector>
#include <map>
#include <functional>
#include <limits>


class KinFit
{  
public:
  KinFit() {}
  ~KinFit() {}
  
  enum class VarProperty_t {
    Gaussian,
    Poissonian,
    LogNormal,
    SquareRoot
  };
  
  void AddMeasuredVariable(const std::string& name, const double value = 0, const double sigma = 0, 
                           const VarProperty_t property = VarProperty_t::Gaussian,
                           const double lowerLimit = std::numeric_limits<double>::quiet_NaN(),
                           const double upperLimit = std::numeric_limits<double>::quiet_NaN()
      );
  void AddUnmeasuredVariable(const std::string& name, const double value = 0);
  void AddFixedVariable(const std::string& name, const double value = 0);
  
  // constraints are functions of the input variables. Should return 0.0 if fulfilled.
  typedef std::function<double (const std::map<std::string, double>&)> constraint_t;
  void AddConstraint(const std::string& name,  constraint_t constraint);
  
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
  std::map<std::string, double> variables; // values of map represent X (works since map is ordered)
  std::vector<double> covariances;         // represents the symmetric covariance matrix
  std::vector<constraint_t> constraints;   // the constraints as C++11 lambdas
  bool initialized;
};

#endif // KINFIT_H
