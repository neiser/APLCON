#ifndef _APLCON_APLCON_HPP
#define _APLCON_APLCON_HPP 1

// detail code is in namespace APLCON_ (note the underscore)
#include "detail/APLCON_hpp.hpp"

#include <algorithm>
#include <functional>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

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

  // the order must correspond to APLCON's status number (see APLCON README)
  enum class Result_Status_t : int {
    Success,
    NoConvergence,
    TooManyIterations,
    UnphysicalValues,
    NegativeDoF,
    OutOfMemory,
    _Unknown // default in Result_t, also used to count items
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

  struct Result_Constraint_t {
    std::string Name;
    size_t Number;    // how many scalar constraints are represented
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
    std::vector<Result_Constraint_t> Constraints;
    const static Result_t Default;
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

  /**
   * @brief DoFit main routine
   * @return the result of the fit, including much additional information
   */
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
                           const double value,
                           const double sigma,
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
                             const double value,
                             const Variable_Settings_t &settings = Variable_Settings_t::Default);
  /**
   * @brief AddFixedVariable
   * @param name unique label for variable
   * @param value initial value for variable
   * @param sigma sqrt of diagonal entry in covariance matrix
   * @param distribution optional type of distribution
   */
  void AddFixedVariable(const std::string& name,
                        const double value,
                        const double sigma,
                        const Distribution_t& distribution = Distribution_t::Gaussian
      );





  void LinkVariable(const std::string& name,
                    const std::vector<double*>& values,
                    const std::vector<double*>& sigmas,
                    const std::vector<Variable_Settings_t>& settings = {}
      );
  void LinkVariable(const std::string& name,
                    const std::vector<double*>& values,
                    const std::vector<double>&  sigmas,
                    const std::vector<Variable_Settings_t>& settings = {}
      );

  void SetCovariance(const std::string& var1, const std::string& var2,
                     const double covariance);
  void SetCovariance(const std::string& var1, const std::string& var2,
                     const std::vector<double>& covariance);
  void SetCovariance(const std::string& var1, const std::string& var2,
                     const std::vector<double*>& covariance);


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
    CheckMapKey("Constraint", name, constraints);

    // define shortcut, but need "typedef", not "using" for older gcc versions...
    typedef APLCON_::function_traits<Functor> trait;

    // non functors are kind of hard to bind later in bind_constraint,
    // so we forbid this here
    static_assert(trait::is_functor, "Only functors are supported as constraints. Wrap and/or bind it if you want to pass such things.");

    using r_type = typename trait::return_type;

    // compile-time check if the Functor returns the proper type
    constexpr bool returns_double = std::is_same<r_type, double >::value;
    constexpr bool returns_vector = std::is_same<r_type, std::vector<double> >::value;
    static_assert(returns_double || returns_vector, "Constraint function does not return double or vector<double>.");

    // compile-time check if the Function wants only double's, or only vector of double's
    // both bool's can never be true at the same time,
    // so we require an exclusive or
    constexpr bool wants_double = trait::template all_args<double>::value;
    constexpr bool wants_vector = trait::template all_args< std::vector<double> >::value;
    static_assert(wants_double ^ wants_vector, "Constraint function does not solely want double or solely vector<double>.");

    // runtime check if given variable number matches to Functor
    constexpr size_t n = trait::arity; // number of arguments in Functor
    if(varnames.size() != n) {
      std::stringstream msg;
      msg << "Constraint '" << name << "': Function argument number (" << n <<
             ") does not match the number of provided varnames (" << varnames.size() << ")";
      throw std::logic_error(msg.str());
    }

    // the flag wants_double and returns_double select the corresponding bind_constraint
    // implementation
    const auto& bound = bind_constraint<returns_double>
        (std::enable_if<wants_double>(),
         constraint, APLCON_::build_indices<n>{});

    constraints[name] = {varnames, bound, wants_double, 0};
    initialized = false;
  }


  // some printout formatting stuff
  // used in overloaded << operators
  // may be changed to tune print formatting
  struct PrintFormatting {
    static std::string Indent;
    static std::string Marker;
    static int Width;
  };

private:

  struct variable_t {
    std::vector<double*> Values;
    std::vector<double*> Sigmas;
    std::vector<Variable_Settings_t> Settings;
    // storage for non-linked variables
    // see AddVariable/LinkVariable for usage
    std::vector<double> StoredValues;
    std::vector<double> StoredSigmas;
    // offsets/indices only valid after Init()
    size_t XOffset;
    std::vector<size_t> V_ij; // for sigmas, see also covariance_t
  };

  struct constraint_t {
    std::vector<std::string> VariableNames;
    std::function< std::vector<double> (const std::vector< std::vector<const double*> >&)> Function;
    bool WantsDouble; // true if Function takes single double as all arguments (set by AddConstraint)
    size_t Number;    // number of represented scalar constraints, set by Init
  };

  // since a variable can represent multiple values
  // we need to store a little symmetric submatrix as a vector (as V) here
  // since we can link values via pointers, it's designed similar to Variable_t
  struct covariance_t {
    std::vector<double*> Values;
    std::vector<double>  StoredValues;
    std::vector<size_t>  V_ij; // indices in V where Values should be copied
  };

  // values with starting values (works since map is ordered)
  typedef std::map<std::string, variable_t> variables_t;
  variables_t variables;
  int nVariables; // number of simple variables
  // off-diagonal covariances addressed by pairs of variable names
  typedef std::map< std::pair<std::string, std::string>, covariance_t > covariances_t;
  covariances_t covariances;
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
  // shortcuts for double limits (used in default values for methods above)
  const static double NaN;

  // private methods
  void Init();
  void AddVariable(const std::string& name, const double value, const double sigma,
                   const APLCON::Variable_Settings_t& settings);

  variables_t::iterator LinkVariable(
      const std::string& name,
      const std::vector<double*>& values,
      const size_t sigma_size,
      const std::vector<Variable_Settings_t>& settings
      );

  covariances_t::iterator MakeCovarianceEntry(
      const std::string& var1, const std::string& var2
      );

  template<typename T>
  void CheckMapKey(const std::string& tag, const std::string& name,
                   std::map<std::string, T> c) {
    if(name.empty()) {
      throw std::logic_error(tag+" name empty");
    }
    if(c.find(name) != c.end()) {
      throw std::logic_error(tag+" with name '"+name+"' already added");
    }
  }

  // define the two different constraint binding functions
  // which are selected on compile-time via their first two arguments

  // the basic idea is to "vectorize" the given constraint function f to fv
  // by defining a lambda fv which is std::bind'ed to the original f
  // then fv can be called on vectors containing pointers to the values
  // on which the constraint should be evaluated
  // see DoFit/Init methods how those arguments for the returned function are constructed

  // is it complicated by the fact that f may return scalar/vector and may want scalar/vector
  // that's why bind_constraint has two dummy arguments which select the correct binding
  // depending on the compile-time analysis of f in AddConstraint. This must be templated because
  // otherwise the compiler evaluates the wrong f call

  using constraint_function_t = std::function< std::vector<double> (const std::vector< std::vector<const double*> >&)>;

  template <bool R, typename F, size_t... I>
  constraint_function_t
  bind_constraint(std::enable_if<true>, // wants_double
                  const F& f, APLCON_::indices<I...>) const {
    auto f_wrap = [] (const F& f, const std::vector< std::vector<const double*> >& x) -> std::vector<double> {
      // dereference the single element inside the inner vector
      // return vector with single element
      return APLCON_::vectorize_if<R>::get(f(*(x[I][0])...));
    };
    return std::move(std::bind(f_wrap, f, std::placeholders::_1));
  }

  template <bool R, typename F, size_t... I>
  constraint_function_t
  bind_constraint(std::enable_if<false>, // wants vector
                  const F& f, APLCON_::indices<I...>) const {
    auto f_wrap = [] (const F& f, const std::vector< std::vector<const double*> >& x) -> std::vector<double> {
      // this might be a little bit inefficient,
      // since we need to allocate the space for the dereferenced double values
      std::vector< std::vector<double> > x_(sizeof...(I));
      for(size_t i=0;i<sizeof...(I);i++) {
        x_[i].resize(x[i].size());
        std::transform(x[i].begin(), x[i].end(),
                       x_[i].begin(),
                       [] (const double* v) { return *v; }
        );
      }
      return APLCON_::vectorize_if<R>::get(f(std::move(x_[I])...));
    };
    return std::move(std::bind(f_wrap, f, std::placeholders::_1));
  }

};

std::ostream& operator<< (std::ostream&, const APLCON::Limit_t&);
std::ostream& operator<< (std::ostream&, const APLCON::Distribution_t&);
std::ostream& operator<< (std::ostream&, const APLCON::Result_Status_t&);
std::ostream& operator<< (std::ostream&, const APLCON::Result_t&);

#endif // APLCON_HPP
