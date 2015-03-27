#ifndef APLCON_HPP
#define APLCON_HPP

#include <vector>
#include <map>
#include <functional>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <typeinfo>
#include <iostream>
#include <sstream>



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
  };

  struct Variable_Settings_t {
    Distribution_t Distribution;
    Limit_t Limit;
    double StepSize;
    const static Variable_Settings_t Default;
  };

  struct Variable_t {
    std::vector<double*> Values;
    std::vector<double> Sigmas;
    std::vector<Variable_Settings_t> Settings;
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

  // the usual constructor
  APLCON(const std::string& _name,
         const Fit_Settings_t& _fit_settings = Fit_Settings_t::Default) :
    instance_name(_name),
    initialized(false),
    instance_id(++instance_counter),
    fit_settings(_fit_settings)
  {}

  // copy the instance,
  // but with new name and possibly new settings
  APLCON(const APLCON& _old,
         const std::string& _name,
         const Fit_Settings_t& _fit_settings = Fit_Settings_t::Default)
    : APLCON(_old)
  {
    instance_name = _name;
    fit_settings  = _fit_settings;
  }


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
                           const double value = NaN,
                           const double sigma = NaN,
                           const Variable_Settings_t &settings = Variable_Settings_t::Default);
  /**
   * @brief AddUnmeasuredVariable
   * @param name unique label for variable
   * @param value initial value for variable
   * @param lowerLimit lower limit of the variable's value
   * @param upperLimit upper limit of the variable's value
   * @param stepSize step size for numerical derivation
   */
  void AddUnmeasuredVariable(const std::string& name,
                             const double value = NaN,
                             const Variable_Settings_t &settings = Variable_Settings_t::Default);
  /**
   * @brief AddFixedVariable
   * @param name unique label for variable
   * @param value initial value for variable
   * @param sigma sqrt of diagonal entry in covariance matrix
   * @param distribution optional type of distribution
   */
  void AddFixedVariable(const std::string& name,
                        const double value = NaN,
                        const double sigma = NaN,
                        const Distribution_t& distribution = Distribution_t::Gaussian
      );


  void SetCovariance(const std::string& var1, const std::string& var2, const double cov);


  void LinkVariable(const std::string& name,
                    const std::vector<double*>& values,
                    const std::vector<double>& sigmas,
                    const std::vector<Variable_Settings_t>& settings = {}
      );


  /**
   * @brief AddConstraint
   * @param name unique label for the constraint
   * @param referred variable names the constraint should act on
   * @param constraint lambda function taking varnames size double arguments, and return double. Should vanish if fulfilled.
   */
  template<typename Functor>
  void AddConstraint(const std::string& name,
                     const std::vector<std::string>& varnames,
                     const Functor& constraint)
  {
    CheckMapKey("Linked constraint", name, constraints);

    using trait = function_traits<Functor>;

    static_assert(trait::is_functor, "Only functors are supported as constraints, wrap and bind it if you want to pass such things.");

    constexpr size_t n = trait::arity; // number of arguments in Functor
    using r_type = typename trait::return_type;

    constexpr bool returns_double = std::is_same<r_type, double >::value;
    constexpr bool returns_vector = std::is_same<r_type, std::vector<double> >::value;
    static_assert(returns_double || returns_vector, "Constraint function does not return double or vector<double>.");

    if(varnames.size() != n) {
      std::stringstream msg;
      msg << "Constraint function argument number (" << n <<
             ") does not match the number of provided varnames (" << varnames.size() << ")";
      throw std::logic_error(msg.str());
    }

    const auto& b = bind_constraint(std::enable_if<returns_double>(),
                                    constraint, build_indices<n>{});

    constraints[name] = {varnames, b};
    initialized = false;
  }


  // some printout formatting stuff
  // used in overloaded << operators
  struct PrintFormatting {
    const static std::string Indent;
    const static std::string Marker;
    const static int Width;
  };


private:

  enum class ConstraintType {
    Simple, Complex
  };

  struct constraint_t {
    std::vector<std::string> VariableNames;
    std::function< std::vector<double> (const std::vector< std::vector<const double*> >&)> Function;
  };


  // values with starting values (works since map is ordered)
  // TODO: merge variables and linked variables
  std::map<std::string, Variable_t> variables;
  int nVariables; // number of simple variables
  // off-diagonal covariances addressed by pairs of variable names
  std::map< std::pair<std::string, std::string>, double > covariances;
  // the constraints
  // a constraint has a list of variable names and
  // a corresponding "vectorized" function evaluated on pointers to double
  std::map<std::string, constraint_t> constraints;
  int nConstraints; // number of double-valued equations, finally determined in Init()

  // storage vectors for APLCON (only usable after Init() call!)
  // X values, V covariances, F constraints
  // and some helper variables
  std::vector<double> X, V, F, V_before;
  std::vector< std::function<std::vector<double>()> > F_func;
  std::map<std::string, size_t> X_s2i; // from varname to index in X

  // since APLCON is stateful, multiple instances of this class
  // need to init APLCON again after switching between them
  // However, when always the same instance is run, we don't need
  // to init APLCON
  std::string instance_name;
  bool initialized;
  static int instance_counter; // global instance counter (never decremented)
  static int instance_lastfit; // save last instance id
  int instance_id;

  // global APLCON settings
  Fit_Settings_t fit_settings;

  void Init();
  void AddVariable(const std::string& name, const double value, const double sigma,
                   const APLCON::Variable_Settings_t& settings);
  template<typename T>
  void CheckMapKey(const std::string& tag, const std::string& name,
                   std::map<std::string, T> c);

  // shortcuts for double limits (used in default values for methods above)
  const static double NaN;

  // some extra stuff for having a nice constraint interface

  // get some function traits like return type and number of arguments
  // based on https://functionalcpp.wordpress.com/2013/08/05/function-traits/

  // functor
  template<class F>
  struct function_traits
  {
  private:
    using call_type = function_traits<decltype(&F::operator())>; // how to prevent the compiler message here for non-callable types?
  public:
    using return_type = typename call_type::return_type;

    static constexpr std::size_t arity = call_type::arity - 1;
    static constexpr bool is_functor = true;

    template <std::size_t N>
    struct argument
    {
      static_assert(N < arity, "error: invalid parameter index.");
      using type = typename call_type::template argument<N+1>::type;
    };

  };

  template<class F>
  struct function_traits<F&> : public function_traits<F> {};

  template<class F>
  struct function_traits<F&&> : public function_traits<F> {};

  // anything else what is callable
  template<typename R, typename... Args>
  struct function_traits<R(Args...)>
  {
    using return_type = R;

    static constexpr std::size_t arity = sizeof...(Args);
    static constexpr bool is_functor = false;

    template <std::size_t N>
    struct argument
    {
      static_assert(N < arity, "error: invalid parameter index.");
      using type = typename std::tuple_element<N,std::tuple<Args...>>::type;
    };

  };

  // function pointer
  template<typename R, typename... Args>
  struct function_traits<R(*)(Args...)> : public function_traits<R(Args...)> {};

  // member function pointer
  template<class C, class R, class... Args>
  struct function_traits<R(C::*)(Args...)> : public function_traits<R(C&,Args...)> {};

  // const member function pointer
  template<class C, class R, class... Args>
  struct function_traits<R(C::*)(Args...) const> : public function_traits<R(C&,Args...)> {};

  // member object pointer
  template<class C, class R>
  struct function_traits<R(C::*)> : public function_traits<R(C&)> {};


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


  // define the two different constraint binding functions
  // which are selected on compile-time via their first argument

  template <typename F, size_t... I>
  std::function< std::vector<double> (const std::vector< std::vector<const double*> >&)>
  bind_constraint(std::enable_if<true>, const F& f, indices<I...>) const {
    // similar to bind_constraint,
    // but f takes now some vector's of doubles (instead of just one double)
    auto fv = [] (const F& f, const std::vector< std::vector<const double*> >& x) -> std::vector<double> {
      // "vectorize" the given constraint function f to fv
      // by defining a lambda fv which is bound to the original f
      // then fv can be called on vectors containing pointers to the values
      // on which the constraint should be evaluated
      // see DoFit/Init methods how those arguments for the returned function are constructed
      std::vector<const double*> x_(sizeof...(I));
      for(size_t i=0;i<sizeof...(I);i++) {
        // size==0 should not appear due to construction in Init
        x_.push_back(x[i][0]);
      }
      // we can simplify now the f call by dereferencing here
      return {f(*(x_[I])...)};
    };
    return std::bind(fv, f, std::placeholders::_1);
  }

  template <typename F, size_t... I>
  std::function< std::vector<double> (const std::vector< std::vector<const double*> >&)>
  bind_constraint(std::enable_if<false>, const F& f, indices<I...>) const {
    // similar to bind_constraint,
    // but f takes now some vector's of doubles (instead of just one double)
    auto fv = [] (const F& f, const std::vector< std::vector<const double*> >& x) {
      // TODO: Make vector<double> instead of pointers...?
      return f(move(x[I])...);
    };

    return std::bind(fv, f, std::placeholders::_1);
  }


};

// templated methods must be implemented in header file

//template<typename T>
//void APLCON::AddConstraint(const std::string& name,
//                   const std::vector<std::string>& varnames,
//                   const T& constraint)
//{
//  CheckMapKey("Constraint", name, constraints);
//  auto f = make_function(constraint);
//  const size_t n = count_arg<decltype(f)>::value;
//  if(varnames.size() != n) {
//    throw std::logic_error("Constraint function argument number does not match the number of varnames.");
//  }
//  constraints[name] = {varnames, bind_constraint(constraint, build_indices<n> {})};
//  initialized = false;
//}

template<typename T>
void APLCON::CheckMapKey(const std::string& tag, const std::string& name, std::map<std::string, T> c) {
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
