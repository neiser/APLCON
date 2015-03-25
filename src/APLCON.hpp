#ifndef APLCON_HPP
#define APLCON_HPP

#include <vector>
#include <map>
#include <functional>
#include <limits>
#include <assert.h>
#include <type_traits>
#include <stdexcept>

/**
 * @brief The APLCON class
 * Provides a C++11'ish wrapper around
 * V.Blobel's FORTRAN APLCON constrained least squares fitter
 * see http://www.desy.de/~blobel/wwwcondl.html for details of the original FORTRAN code
 */
class APLCON
{
public:

  struct Fit_Settings_t {
    int DebugLevel;
    int MaxIterations;
    double ConstraintAccuracy;
    double MeasuredStepSizeFactor;
    double UnmeasuredStepSizeFactor;
    double MinimalStepSizeFactor;
    const static Fit_Settings_t Default;
  };

  enum class Distribution_t {
    Gaussian,
    Poissonian,
    LogNormal,
    SquareRoot
  };


  struct Limit_t {
    double Low;
    double High;
    const static Limit_t NoLimit;
  };

  struct Variable_Settings_t {
    Distribution_t Distribution;
    Limit_t Limit;
    double StepSize;
  };

  struct Variable_t {
    double Value;
    double Sigma;
    Variable_Settings_t Settings;
  };

  enum class Result_Status_t {
    Success,
    NoConvergence,
    TooManyIterations,
    UnphysicalValues,
    NegativeDoF,
    OutOfMemory
  };

  template<typename T>
  struct Result_BeforeAfter_t {
    T Before;
    T After;
  };


  struct Result_Variable_t {
    std::string Name;
    Result_BeforeAfter_t<double> Value;
    Result_BeforeAfter_t<double> Sigma;
    Result_BeforeAfter_t< std::vector<double> > Covariances;
    double Pull;
    Variable_Settings_t Settings;
  };

  struct Result_t {
    std::string Name;
    Result_Status_t Status;
    double ChiSquare;
    int NDoF;
    double Probability;
    int NIterations;
    int NFunctionCalls;
    std::vector<Result_Variable_t> Variables;
    std::vector<std::string> Constraints;
  };


  APLCON(const std::string& _name = "",
         const Fit_Settings_t& _fit_settings = Fit_Settings_t::Default) :
    instance_name(_name),
    initialized(false),
    instance_id(++instance_counter),
    fit_settings(_fit_settings)
  {}

  Result_t DoFit();

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
                           const Limit_t limit = Limit_t::NoLimit,
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
                             const Limit_t limit = Limit_t::NoLimit,
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
   * @param name unique label for the constraint
   * @param referred variable names the constraint should act on
   * @param constraint lambda function taking varnames size double arguments, and return double. Should vanish if fulfilled.
   */
  template<typename T>
  void AddConstraint(const std::string& name,
                     const std::vector<std::string>& varnames,
                     T constraint);

  void AddCovariance(const std::string& var1, const std::string& var2, const double cov);

  // some printout formatting stuff
  // used in overloaded << operators
  struct PrintFormatting {
    const static std::string Indent;
    const static std::string Marker;
    const static int Width;
  };

private:
  struct constraint_t {
    std::vector<std::string> VariableNames;
    std::function<double(const std::vector<const double*>&)> Function;
  };

  // values with starting values (works since map is ordered)
  std::map<std::string, Variable_t> variables;
  // off-diagonal covariances addressed by pairs of variable names
  std::map< std::pair<std::string, std::string>, double > covariances;
  // the constraints
  // a constraint has a list of variable names and
  // a corresponding "vectorized" function evaluated on pointers to double
  std::map<std::string, constraint_t> constraints;

  // storage vectors for APLCON (only usable after Init() call!)
  // X values, V covariances, F constraints
  // and some helper variables
  std::vector<double> X, V, F, V_before;
  std::vector< std::function<double()> > F_func;
  std::map<std::string, size_t> X_s2i; // from varname to index in X

  // since APLCON is stateful, multiple instances of this class
  // need to init APLCON again after switching between them
  // However, when always the same instance is run, we don't need
  // to init APLCON
  const std::string instance_name;
  bool initialized;
  static int instance_counter; // global instance counter (never decremented)
  static int instance_lastfit; // save last instance id
  int instance_id;

  // global APLCON settings
  Fit_Settings_t fit_settings;

  void Init();
  void AddVariable(const std::string& name, const double value, const double sigma,
                   const Distribution_t distribution,
                   const Limit_t limit,
                   const double stepSize
                   );
  template<typename T>
  void TestName(const std::string& tag, const std::string& name,
                std::map<std::string, T> c);

  // some extra stuff for having a nice constraint interface

  // first it seems pretty complicated to figure out how many arguments a
  // given lambda has (std::function is easy though)
  // based on http://stackoverflow.com/questions/20722918/how-to-make-c11-functions-taking-function-parameters-accept-lambdas-automati/
  // and http://stackoverflow.com/questions/9044866/how-to-get-the-number-of-arguments-of-stdfunction

  template <typename T>
  struct function_traits
     : public function_traits<decltype(&T::operator())>
  {};

  template <typename ClassType, typename ReturnType, typename... Args>
  struct function_traits<ReturnType(ClassType::*)(Args...) const> {
     typedef std::function<ReturnType (Args...)> f_type;
  };

  template <typename L>
  typename function_traits<L>::f_type make_function(L l) {
    return (typename function_traits<L>::f_type)(l);
  }

  // this works only with std::function due to the necessary return type R (I guess)
  template<typename T>
  struct count_arg;

  template<typename R, typename... Args>
  struct count_arg<std::function<R(Args...)>> {
      static const size_t value = sizeof...(Args);
  };

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

// templated methods must be implemented in header file

template<typename T>
void APLCON::AddConstraint(const std::string& name,
                   const std::vector<std::string>& varnames,
                   T constraint)
{
  TestName("Constraint", name, constraints);
  auto f = make_function(constraint);
  const size_t n = count_arg<decltype(f)>::value;
  assert(varnames.size() == n);
  constraints[name] = {varnames, bind_constraint(constraint, build_indices<n> {})};
  initialized = false;
}

template<typename T>
void APLCON::TestName(const std::string& tag, const std::string& name, std::map<std::string, T> c) {
  if(name.empty()) {
    throw std::logic_error(tag+" name empty");
  }
  if(c.find(name) != c.end()) {
    throw std::logic_error(tag+" with name '"+name+"' already added");
  }
}

std::ostream& operator<< (std::ostream&, const APLCON::Limit_t&);
std::ostream& operator<< (std::ostream&, const APLCON::Distribution_t&);
std::ostream& operator<< (std::ostream&, const APLCON::Result_Status_t&);
std::ostream& operator<< (std::ostream&, const APLCON::Result_t&);


#endif // APLCON_HPP
