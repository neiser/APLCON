#ifndef KINFIT_H
#define KINFIT_H

#include <vector>
#include <map>
#include <functional>
#include <limits>
#include <assert.h>
#include <type_traits>

/**
 * @brief The KinFit class
 * Provides a C++11'ish wrapper around
 * V.Blobel's FORTRAN APLCON constrained least squares fitter
 * see http://www.desy.de/~blobel/wwwcondl.html for details
 */
class KinFit
{
public:
  KinFit() : initialized(false), constraint_added(false) {}
  ~KinFit() {}

  enum class Distribution_t {
    Gaussian,
    Poissonian,
    LogNormal,
    SquareRoot
  };

  using map_t = std::map<std::string, double>;

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

  /**
   * @brief AddConstraint
   *
   * @param name
   * @param f
   */
  template<typename varnames_t, size_t n_varnames, typename constraint_function_t>
  void AddConstraint(const std::string& name, varnames_t const(&varnames)[n_varnames], constraint_function_t constraint)
  {
    // convert the vars into vector string
    // we need to have this weird C-style array to deduce the size N of parameters
    static_assert(std::is_same<varnames_t, const char*>::value, "AddConstraint only works with constant strings as variable names");
    const std::vector<std::string> v(std::begin(varnames), std::end(varnames));
    constraints[name] = {v, bind_constraint(constraint, build_indices<n_varnames>{})};
  }

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
  // values with starting values (works since map is ordered)
  std::map<std::string, double> variables;
  // track the type of distributions
  std::vector<Distribution_t> distributions;
  // track the limits (NaN if unset)
  std::vector< std::pair<double, double> > limits;
  // represents the symmetric covariance matrix
  std::vector<double> covariances;
  // the constraints
  // a constraint has a list of variable names and
  // a corresponding "vectorized" function evaluated on pointers to double
  struct constraint_t {
    std::vector<std::string> VariableNames;
    std::function<double(const std::vector<const double*>&)> Function;
  };
  std::map<std::string, constraint_t> constraints;
  // step sizes for numerical evaluation (zero if fixed, NaN if APLCON default)
  std::vector<double> stepSizes;

  bool initialized;
  bool constraint_added;

  void Init();
  void AddVariable(const std::string& name, const double value, const double sigma,
                   const Distribution_t distribution,
                   const double lowerLimit, const double upperLimit,
                   const double stepSize
                   );

  // some extra stuff for having a nice constraint interface

  // this little template fun is called "pack of indices"
  // it enables the nice definition of constraints via AddConstraint(...) method
  // see http://stackoverflow.com/questions/11044504/any-solution-to-unpack-a-vector-to-function-arguments-in-c
  // and http://loungecpp.wikidot.com/tips-and-tricks%3aindices
  template <std::size_t... Is>
  struct indices {};

  template <std::size_t N, std::size_t... Is>
  struct build_indices : build_indices<N-1, N-1, Is...> {};

  template <std::size_t... Is>
  struct build_indices<0, Is...> : indices<Is...> {};

  template <typename FuncType, size_t... I>
  std::function<double(const std::vector<const double*>&)> bind_constraint(const FuncType& f, indices<I...>) const {
    // "vectorize" the given constraint function f to fv
    // by defining a lambda fv which is bound to the original f
    // then fv can be called on vectors containing pointers to the values
    // on which the constraint should be evaluated
    // see DoFit how this is actually done
    auto fv = [] (const FuncType& f, const std::vector<const double*>& x) {
      assert(x.size() == sizeof...(I));
      return f(*x[I]...);
    };
    return std::bind(fv, f, std::placeholders::_1);
  }

};

#endif // KINFIT_H
