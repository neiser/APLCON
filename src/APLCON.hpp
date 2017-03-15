#ifndef _APLCON_APLCON_HPP
#define _APLCON_APLCON_HPP 1

// detail code is in namespace APLCON_ (note the underscore)
#include "detail/APLCON_hpp.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

/**
 * @brief The APLCON C++11 wrapper class
 *
 * It provides a C++11'ish wrapper around
 * V.Blobel's constrained least squares fitter written in Fortran.
 * See http://www.desy.de/~blobel/wwwcondl.html for details of the original FORTRAN code.
 */

class APLCON
{
public:

   /**
   * @brief The Fit_Settings_t struct. See APLCON itself for details.
   */
  struct Fit_Settings_t {
    int DebugLevel;
    int MaxIterations;
    double ConstraintAccuracy;
    double Chi2Accuracy;
    double MeasuredStepSizeFactor;
    double UnmeasuredStepSizeFactor;
    double MinimalStepSizeFactor;
    bool   SkipCovariancesInResult;
    const static Fit_Settings_t Default;
  };

  /**
   * @brief The Distribution_t enum
   * @see Variable_Settings_t
   */
  enum class Distribution_t {
    Gaussian, /**< Gaussian distributed variable (default) */
    Poissonian, /**< Poissonian distributed variable */
    LogNormal, /**< Ratios are lognormal distributed */
    SquareRoot /**< SquareRoot transformation */
  };

  /**
   * @brief The Limit_t struct defines upper and lower limits
   */
  struct Limit_t {
    double Low;
    double High;
  };

  /**
   * @brief The Variable_Settings_t struct contains settings per variable
   */
  struct Variable_Settings_t {
    Distribution_t Distribution;
    Limit_t Limit;
    double StepSize;
    const static Variable_Settings_t Default;
  };

  /**
   * @brief The Result_Status_t enum encodes APLCON's status after fit
   * @note the order must correspond to APLCON's status number (see APLCON README)
   */
  enum class Result_Status_t : int {
    Success, /**< Fit was successful, result is meaningful */
    NoConvergence, /**< No convergence reached */
    TooManyIterations, /**< Too many iterations needed with no convergence */
    UnphysicalValues, /**< Unphysical values encountered during fit */
    NegativeDoF, /**< Negative degrees of freedom */
    OutOfMemory, /**< Not sufficient memory for fit */
    _Unknown // default in Result_t, also used to count items
  };

  /**
   * @brief The Result_BeforeAfter_t struct wraps information before and after the fit
   * @see Result_Variable_t
   */
  template<typename T>
  struct Result_BeforeAfter_t {
    T Before;
    T After;
  };

  /**
   * @brief The Result_Variable_t struct contains the fit result per variable
   */
  struct Result_Variable_t {
    std::string PristineName; /**< pristine name, i.e. without appended "[i]" */
    size_t Dimension; /**< dimension of vector variable, =1 for scalar variable */
    size_t Index;     /**< index inside vector variable, =0 for scalar variable */
    Result_BeforeAfter_t<double> Value;
    Result_BeforeAfter_t<double> Sigma;
    Result_BeforeAfter_t< std::map<std::string, double> > Covariances;
    double Pull;
    Variable_Settings_t Settings;
  };

  /**
   * @brief helper method to calculate correlations from covariances
   * @note used by the ostream<< operators, and might be useful for some users
   * @param variables from Result_t
   * @return correlations between stringified variables, before and after fit
   */
  static std::map< std::string, Result_BeforeAfter_t< std::map<std::string, double> > >
  CalculateCorrelations(const std::map<std::string, Result_Variable_t>& variables);

  /**
   * @brief The Result_Constraint_t struct contains constraint information
   */
  struct Result_Constraint_t {
    size_t Dimension;    // how many scalar constraints are represented by it
  };

  /**
   * @brief The Result_t struct contains
   * after the fit all information about it.
   * Variables are always referenced by their string representation,
   * appended with [i] if they are non-scalar.
   */
  struct Result_t {
    std::string Name;
    Result_Status_t Status;
    double ChiSquare;
    int NDoF;
    double Probability;
    int NIterations;
    int NFunctionCalls;
    std::map<std::string, Result_Variable_t>   Variables;
    std::map<std::string, Result_Constraint_t> Constraints;
    int NScalarConstraints;
    const static Result_t Default;
  };

  /**
   * @brief Create new APLCON instance with a name, and optional fit settings
   * @param _name
   * @param _fit_settings
   */
  APLCON(const std::string& _name,
         const Fit_Settings_t& _fit_settings = Fit_Settings_t::Default) :
    instance_name(_name),
    initialized(false),
    instance_id(++instance_counter),
    fit_settings(_fit_settings) {}

  /**
   * @brief Copy instance from an existing one with new name and new settings
   * @param _old instance to be copied from
   * @param _name name of new instance
   * @param _fit_settings new fit settings
   */
  APLCON(const APLCON& _old,
         const std::string& _name,
         const Fit_Settings_t& _fit_settings)
    : APLCON(_old)
  {
    instance_name = _name;
    fit_settings  = _fit_settings;
  }

  /**
   * @brief Copy instance from an existing one with new name
   * @param _old instance to be copied from
   * @param _name name of new instance
   */
  APLCON(const APLCON& _old,
         const std::string& _name)
    : APLCON(_old, _name, _old.fit_settings) {}

  /**
   * @brief Get the Name of this APLCON instance
   * @return the name
   */
  std::string GetName() const {
      return instance_name;
  }

  /**
   * @brief Obtain current fitter settings
   * @return settings to be obtained
   */
  const Fit_Settings_t& GetSettings() const {
    return fit_settings;
  }

  /**
   * @brief Set given fitter settings
   * @param _new_settings new settings struct
   * @note the fitter is un-initialized after the call of this method
   */
  void SetSettings(const Fit_Settings_t& _new_settings) {
    initialized = false;
    fit_settings = _new_settings;
  }

  /**
   * @brief Obtain variable names
   * @return vector of build variable names which have been added so far
   */
  std::vector<std::string> VariableNames() const;

  /**
   * @brief Main routine of the fitter
   * @return the result of the fit, including all information
   */
  Result_t DoFit();

  /**
   * @brief Add measured variable to fitter with given sigma, internally stored
   * @see LinkVariable for linking externally stored values
   * @param name unique label for variable
   * @param value initial value for variable
   * @param sigma sqrt of diagonal entry in covariance matrix
   * @param settings additional settings like distribution or limits
   */
  void AddMeasuredVariable(const std::string& name,
                           double value,
                           double sigma,
                           const Variable_Settings_t &settings = Variable_Settings_t::Default);
  /**
   * @brief Add unmeasured value to fitter which has vanishing sigma
   * @param name unique label for variable
   * @param value initial value for variable
   * @param settings additional settings like distribution or limits
   */
  void AddUnmeasuredVariable(const std::string& name,
                             double value = 0,
                             const Variable_Settings_t &settings = Variable_Settings_t::Default);
  /**
   * @brief Add a fixed variable to fitter
   * @param name unique label for variable
   * @param value initial value for variable
   * @param sigma sqrt of diagonal entry in covariance matrix
   * @param distribution optional type of distribution, default is Gaussian
   */
  void AddFixedVariable(const std::string& name,
                        double value,
                        double sigma,
                        const Distribution_t& distribution = Distribution_t::Gaussian
      );


  /**
   * @brief Link externally stored variable to fitter
   * @param name unique label for variable
   * @param values vector of double pointers, owned by user
   * @param sigmas corresponding sigmas, provide zero for unmeasured value
   * @param settings corresponding optional settings, provide stepsize=0 for fixed
   */
  void LinkVariable(const std::string& name,
                    const std::vector<double*>& values,
                    const std::vector<double*>& sigmas,
                    const std::vector<Variable_Settings_t>& settings = DefaultSettings
      );

  /**
   * @brief LinkVariable with values only, sigmas are internally stored
   * @param name unique label for variable
   * @param values vector of double pointers, owned by user
   * @param sigmas vector of double values, managed by APLCON instance
   * @param settings optional settings
   */
  void LinkVariable(const std::string& name,
                    const std::vector<double*>& values,
                    const std::vector<double>&  sigmas,
                    const std::vector<Variable_Settings_t>& settings = DefaultSettings
      );

  /**
   * @brief LinkVariable with values, sigmas and pulls
   * @param name unique label for variable
   * @param values starting values, overwritten by result after DoFit()
   * @param sigmas pointers to uncertainties, provide zero for unmeasured
   * @param pulls designated output for pulls after DoFit() (given values are ignored/overwritten)
   * @param settings optional settings for each variable
   */
  void LinkVariable(const std::string& name,
                    const std::vector<double*>& values,
                    const std::vector<double*>& sigmas,
                    const std::vector<double*>& pulls,
                    const std::vector<Variable_Settings_t>& settings = DefaultSettings
      );

  /**
   * @brief Set covariance between to variable names
   * @param var1 first variable name
   * @param var2 second variable name
   * @param covariances with correct size, see advanced example. May provide APLCON::NaN to ignore.
   */
  void SetCovariance(const std::string& var1, const std::string& var2,
                     const std::vector<double>& covariances);

  /**
   * @brief Set covariance between two scalar variables
   * @param var1 first variable name
   * @param var2 second variable name
   * @param covariance one value for two scalar variables var1 and var2
   */
  void SetCovariance(const std::string& var1, const std::string& var2,
                     double covariance);

  /**
   * @brief Link covariance value between two variable names
   * @param var1 first variable name
   * @param var2 second variable name
   * @param covariances vector of pointers with correct size (may provide nullptr to ignore)
   */
  void LinkCovariance(const std::string& var1, const std::string& var2,
                      const std::vector<double*>& covariances);


  /**
   * @brief Add named constraint to the fitter depending on given variables
   * @param name unique label for the constraint
   * @param varnames variable names the constraint should act on
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
    constexpr size_t n = trait::arity; // number of arguments in Functor
    constexpr bool wants_double = trait::template all_args<double>::value;
    constexpr bool wants_vector = trait::template all_args< std::vector<double> >::value;
    constexpr bool wants_matrix = n==1 && trait::template all_args< std::vector< std::vector<double> > >::value;
    static_assert(wants_double + wants_vector + wants_matrix == 1,
                  "Constraint function does not either take double's, or vector<double>'s, or single vector<vector<double>> (matrix) as argument(s).");


    // runtime check if given variable number matches to Functor
    if(!wants_matrix && varnames.size() != n) {
      std::stringstream msg;
      msg << "Constraint '" << name << "': Function argument number (" << n <<
             ") does not match the number of provided varnames (" << varnames.size() << ")";
      throw Error(msg.str());
    }

    // the flag wants_double and returns_double select the corresponding bind_constraint
    // implementation
    const auto& bound = bind_constraint<returns_double>
        (std::enable_if<wants_double>(),
         std::enable_if<wants_vector>(),
         constraint, APLCON_::build_indices<n>{});

    constraints[name] = {varnames, bound, wants_double, 0};
    initialized = false;
  }

  // shortcuts for double limits (used in default values for methods above)
  constexpr static double NaN = std::numeric_limits<double>::quiet_NaN(); /**< short cut for NaN value */
  static std::vector<Variable_Settings_t> DefaultSettings; /**< short cut for empty variable settings */

  /**
   * @brief The Error class is thrown if anything is not correct in the fitter setup
   */
  class Error : public std::runtime_error {
  public:
    explicit Error(const std::string& msg) : runtime_error(msg) {}
  };

  /**
   * @brief The PrintFormatting struct tunes the output of ostream<< operators for APLCON::Result_t
   */
  struct PrintFormatting {
    static std::string Indent;
    static std::string Marker;
    static int Width;
  };

private:

  // instances should be copied so easily
  // you're required to give the new instance a name
  // (which is hopefully different to the one before)
  APLCON(const APLCON&) = default;
  APLCON& operator= (const APLCON& other) = delete;

  struct variable_t {
    std::vector<double*> Values;
    std::vector<double*> Sigmas;
    std::vector<double*> Pulls;
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

  // private methods
  void Init();
  void InitAPLCON();
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

  APLCON::variable_t GetVariableByName(const std::string& varname, const std::string& errmsg) {
    const auto& it = variables.find(varname);
    if(it == variables.end()) {
      throw Error(errmsg);
    }
    return it->second;
  }

  template<typename T>
  void CheckMapKey(const std::string& tag, const std::string& name,
                   std::map<std::string, T> c) {
    if(name.empty()) {
      throw Error(tag+" name empty");
    }
    if(c.find(name) != c.end()) {
      throw Error(tag+" with name '"+name+"' already added");
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
  bind_constraint(std::enable_if<true>,  // wants double
                  std::enable_if<false>, // does not want vector
                  const F& f, APLCON_::indices<I...>) const {
    return [f] (const std::vector< std::vector<const double*> >& x) -> std::vector<double> {
      // dereference the single element inside the inner vector
      // return vector with single element
      return APLCON_::vectorize_if<R>::get(f(*(x[I][0])...));
    };
  }

  template <bool R, typename F, size_t... I>
  constraint_function_t
  bind_constraint(std::enable_if<false>, // does not want double
                  std::enable_if<true>,  // wants vector
                  const F& f, APLCON_::indices<I...>) const {
    return [f] (const std::vector< std::vector<const double*> >& x) -> std::vector<double> {
      // this might be a little bit inefficient,
      // since we need to allocate the space for the dereferenced double values
      // but well, the constraints then look easier
      std::vector< std::vector<double> > x_(sizeof...(I));
      for(size_t i=0;i<sizeof...(I);i++) {
        x_[i].resize(x[i].size());
        std::transform(x[i].begin(), x[i].end(),
                       x_[i].begin(),
                       [] (const double* v) { return *v; }
        );
      }
      return APLCON_::vectorize_if<R>::get(f(x_[I]...));
    };
  }

  template <bool R, typename F, size_t... I>
  constraint_function_t
  bind_constraint(std::enable_if<false>, // does not want double
                  std::enable_if<false>, // does not want vector, so wants matrix!
                  const F& f, APLCON_::indices<I...>) const {
    return [f] (const std::vector< std::vector<const double*> >& x) -> std::vector<double> {
      // this might be a little bit inefficient,
      // since we need to allocate the space for the dereferenced double values
      // but well, the constraints then look easier
      std::vector< std::vector<double> > x_(x.size());
      for(size_t i=0;i<x.size();i++) {
        x_[i].resize(x[i].size());
        std::transform(x[i].begin(), x[i].end(),
                       x_[i].begin(),
                       [] (const double* v) { return *v; }
        );
      }
      return APLCON_::vectorize_if<R>::get(f(x_));
    };
  }

};

/** @example src/example/00_verysimple.cc */
/** @example src/example/01_simple.cc */
/** @example src/example/02_linker.cc */
/** @example src/example/03_advanced.cc */

std::ostream& operator<< (std::ostream&, const APLCON::Limit_t&);
std::ostream& operator<< (std::ostream&, const APLCON::Distribution_t&);
std::ostream& operator<< (std::ostream&, const APLCON::Result_Status_t&);
std::ostream& operator<< (std::ostream&, const APLCON::Result_t&);

#endif // APLCON_HPP
