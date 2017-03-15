#include "APLCON.hpp"

// detail code is in namespace APLCON_ (note the underscore)
#include "detail/APLCON_cc.hpp"

// long ostream stuff is in extra header
#include <detail/APLCON_ostream.hpp>

// include the wrapping C functions
extern "C" {
#include "wrapper/APLCON.h"
}

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

using namespace std;

int APLCON::instance_counter = 0;
int APLCON::instance_lastfit = 0;
std::vector<APLCON::Variable_Settings_t> APLCON::DefaultSettings;

const APLCON::Variable_Settings_t APLCON::Variable_Settings_t::Default = {
  APLCON::Distribution_t::Gaussian,
  {
    -numeric_limits<double>::infinity(),
    numeric_limits<double>::infinity()
  },
  APLCON::NaN
};

// transferred to APLCON in Init() method
const APLCON::Fit_Settings_t APLCON::Fit_Settings_t::Default = {
  0,  // no debug printout from APLCON
  -1, // max iterations
  APLCON::NaN, // ConstraintAccuracy
  APLCON::NaN, // Chi2Accuracy
  APLCON::NaN, // MeasuredStepSizeFactor
  APLCON::NaN, // UnmeasuredStepSizeFactor
  APLCON::NaN, // MinimalStepSizeFactor
  false,       // SkipCovariancesInResult
};

// proper default result
const APLCON::Result_t APLCON::Result_t::Default = {
  "", // Name
  APLCON::Result_Status_t::_Unknown,
  APLCON::NaN,
  -1,
  APLCON::NaN,
  -1,
  -1,
  {},
  {},
  -1
};

// simple AddVariable methods

void APLCON::AddMeasuredVariable(const std::string &name, const double value, const double sigma,
                                 const APLCON::Variable_Settings_t& settings)
{
  if(sigma == 0) {
    throw Error("Measured variables need non-zero sigma. By definition, they are unmeasured then.");
  }
  if(settings.StepSize == 0) {
    throw Error("Measured variables need non-zero step size. By definition, they are fixed then.");
  }
  AddVariable(name, value, sigma, settings);
}

void APLCON::AddUnmeasuredVariable(const string &name, const double value,
                                   const APLCON::Variable_Settings_t& settings)
{
  if(settings.StepSize == 0) {
    throw Error("Unmeasured variables need non-zero step size. By definition, they are fixed then.");
  }
  // unmeasured variables have a sigma of 0
  AddVariable(name, value, 0, settings);
}

void APLCON::AddFixedVariable(const string &name, const double value, const double sigma,
                              const Distribution_t &distribution)
{
  if(sigma == 0) {
    throw Error("Fixed variables need non-zero sigma. By definition, they are unmeasured then.");
  }
  // fixed variables have stepSize of 0
  // and limits don't apply (probably?)
  AddVariable(name, value, sigma, {distribution, Variable_Settings_t::Default.Limit, 0});
}

void APLCON::AddVariable(const string &name, const double value, const double sigma,
                         const APLCON::Variable_Settings_t& settings)
{
  // check if variable already exists
  CheckMapKey("Variable", name, variables);

  // only allow finite low and high limits
  if(isfinite(settings.Limit.Low) ^ isfinite(settings.Limit.High)) {
    throw Error("Variable '"+name+"' does not specify High AND Low limit");
  }

  variable_t var;
  var.StoredValues.emplace_back(value);
  var.StoredSigmas.emplace_back(sigma);
  var.Settings.emplace_back(settings);

  variables[name] = var;
  initialized = false;
}

// LinkVariable methods

void APLCON::LinkVariable(const string &name,
                          const vector<double*> &values,
                          const vector<double*> &sigmas,
                          const vector<APLCON::Variable_Settings_t> &settings) {
  // linked sigmas are easy
  auto it = LinkVariable(name, values, sigmas.size(), settings);
  it->second.Sigmas = sigmas;
}

void APLCON::LinkVariable(const string &name,
                          const vector<double*>& values,
                          const vector<double>& sigmas,
                          const vector<APLCON::Variable_Settings_t> &settings) {
  const size_t sigmas_size = sigmas.size()==1 ? values.size() : sigmas.size();
  auto it = LinkVariable(name, values, sigmas_size, settings);
  it->second.StoredSigmas = sigmas;
  if(sigmas.size()==1) {
    it->second.StoredSigmas.resize(sigmas_size, sigmas[0]);
  }
}

void APLCON::LinkVariable(const string& name,
                          const std::vector<double*>& values,
                          const std::vector<double*>& sigmas,
                          const std::vector<double*>& pulls,
                          const std::vector<APLCON::Variable_Settings_t>& settings) {
  if(values.size() != pulls.size()) {
    throw Error("Pulls size does not match number of provided values");
  }
  auto it = LinkVariable(name, values, sigmas.size(), settings);
  it->second.Sigmas = sigmas;
  it->second.Pulls = pulls;
}

APLCON::variables_t::iterator APLCON::LinkVariable(
    const string &name,
    const std::vector<double*>& values,
    const size_t sigmas_size,
    const std::vector<APLCON::Variable_Settings_t> &settings
    ) {

  CheckMapKey("Linked Variable", name, variables);

  const size_t n = values.size();
  if(n==0) {
    throw Error("At least one value should be linked");
  }

  if(sigmas_size != n) {
    throw Error("Sigmas size does not match number of provided values");
  }

  variable_t var;

  // check settings, set to defaults if none given
  if(settings.size() == 0) {
    var.Settings.resize(n, Variable_Settings_t::Default);
  }
  else if(settings.size() == 1) {
    var.Settings.resize(n, settings[0]);
  }
  else if(settings.size() != n) {
    // 0, 1 or it must fit to values...
    throw Error("Settings size does not match number of provided values");
  }
  else {
    var.Settings = settings;
  }

  // set values, that's easy
  var.Values = values;

  auto p = variables.insert(make_pair(name, var));
  return p.first;
}

// SetCovariance methods

// assume var1 and var2 are both scalar
void APLCON::SetCovariance(const std::string& var1,
                           const std::string& var2,
                           const double covariance) {

  if(var1 == var2) {
    throw Error("Covariance variable names must be different");
  }
  auto it = MakeCovarianceEntry(var1, var2);
  auto& values = it->second.StoredValues;
  values.resize(1);
  values.back() = covariance;
}
void APLCON::SetCovariance(const std::string& var1,
                           const std::string& var2,
                           const std::vector<double>& covariances) {
  auto it = MakeCovarianceEntry(var1, var2);
  if(covariances.empty()) {
    throw Error("Empty covariance values given");
  }
  it->second.StoredValues = covariances;
}

void APLCON::LinkCovariance(const std::string& var1,
                           const std::string& var2,
                           const std::vector<double*>& covariances) {
  auto it = MakeCovarianceEntry(var1, var2);
  if(covariances.empty()) {
    throw Error("Empty covariance pointers given");
  }
  it->second.Values = covariances;
}

APLCON::covariances_t::iterator APLCON::MakeCovarianceEntry(
    const string &var1,
    const string &var2)
{
  if(var1.empty() || var2.empty()) {
    throw Error("Covariance variable names cannot be empty strings");
  }
  // we cannot check in general if the variables should have different names
  // since for vector variables, there are still correlations possible
  // finally, this will be checked by Init

  // assume covariance is gonna changed
  // so force Init to execute again
  initialized = false;

  // search the covariances map
  // but note that the pairs are symmetric (since the covariance matrix is)
  // so we search for both possibilities
  // if we find it, we update it...
  const pair<string, string>& p1 = make_pair(var1, var2);
  auto it1 = covariances.find(p1);
  if(it1 != covariances.end()) {
    return it1;
  }
  const pair<string, string>& p2 = make_pair(var2, var1);
  auto it2 = covariances.find(p2);
  if(it2 != covariances.end()) {
    return it2;
  }

  // not found, then add default struct and return
  auto p = covariances.insert(make_pair(p1, covariance_t()));
  return p.first;
}


map< string, APLCON::Result_BeforeAfter_t< map<string, double> > >
APLCON::CalculateCorrelations(const map<string, Result_Variable_t>& variables)
{
  // slow but simple implementation
  // see covariances building in DoFit for a faster implementation

  map< string, APLCON::Result_BeforeAfter_t< map<string, double> > > correlations;

  for(const auto& it_map_i : variables) {

    for(const auto& it_map_j : variables) {

      const string& varname_i = it_map_i.first;
      const Result_Variable_t& var_i = it_map_i.second;
      const string& varname_j = it_map_j.first;
      const Result_Variable_t& var_j = it_map_j.second;

      if(var_i.Covariances.Before.count(varname_i) == 0)
          continue;
      if(var_j.Covariances.Before.count(varname_j) == 0)
          continue;

      const double prod_before =
          var_i.Covariances.Before.at(varname_i) *
          var_j.Covariances.Before.at(varname_j);
      const double prod_after  =
          var_i.Covariances.After.at(varname_i) *
          var_j.Covariances.After.at(varname_j);

      correlations[varname_i].Before[varname_j] =
          var_i.Covariances.Before.at(varname_j) / sqrt(prod_before);

      correlations[varname_i].After[varname_j] =
          var_i.Covariances.After.at(varname_j) / sqrt(prod_after);
    }
  }

  return correlations;
}

vector<string> APLCON::VariableNames() const {
  vector<string> variableNames;
  for(const auto& it_map : variables) {
    const string& name = it_map.first;
    const variable_t& var = it_map.second;
    // figuring out the number of variables here
    // is difficult since Values might not be built by Init yet
    size_t n = var.Values.size()==0 ? var.StoredValues.size() : var.Values.size();
    for(size_t k=0;k<n;k++) {
      const string& varname = APLCON_::BuildVarName(name, n, k);
      variableNames.push_back(varname);
    }
  }
  return variableNames;
}

// Main Fit Routines


APLCON::Result_t APLCON::DoFit()
{
  // ensure that APLCON is properly initialized
  // and the value vectors X, F, V
  Init();

  // the main convergence loop
  int aplcon_ret = -1;
  do {
    // evaluate the constraints F_func and
    // store results in F via iterator F_it
    auto F_it = F.begin();
    for(size_t i=0; i<F_func.size(); ++i) {
      for(const auto& v : F_func[i]()) {
        *F_it = v;
        ++F_it;
      }
    }
    // call APLCON iteration
    c_aplcon_aploop(X.data(), V.data(), F.data(), &aplcon_ret);
  }
  while(aplcon_ret<0);

  // now, after the loop, retrieve the results
  Result_t result = Result_t::Default;

  // make some evil static_cast, but it's way shorter than switch statement
  if(aplcon_ret >= static_cast<int>(Result_Status_t::_Unknown)) {
    throw Error("Unkown return value after APLCON fit");
  }
  result.Status = static_cast<Result_Status_t>(aplcon_ret);

  // now retrieve "everything" from APLCON

  // retrieve some info about the fit (directly copy to struct field if possible)
  float chi2, pval;
  // chndpv and apstat both return the resulting chi2,
  // but the latter returns it with double precision
  c_aplcon_chndpv(&chi2,&result.NDoF,&pval);
  result.Probability = pval;
  c_aplcon_apstat(&result.ChiSquare, &result.NFunctionCalls, &result.NIterations);

  // get the pulls from APLCON
  vector<double> pulls(X.size());
  c_aplcon_appull(pulls.data());

  // copy just the names of the constraints
  for(const auto& it_map : constraints) {
    Result_Constraint_t r_con;
    r_con.Dimension = it_map.second.Number;
    result.Constraints[it_map.first] = r_con;
  }
  result.NScalarConstraints = nConstraints;
  result.Name = instance_name;

  // now we're ready to fill the result.Variables vector
  // we fill it in the same order as the X vector
  // which makes debug output from  APLCON comparable to dumps of this structure

  for(const auto& it_map : variables) {
    const string& name = it_map.first;
    const variable_t& before = it_map.second;

    for(size_t k=0;k<before.Values.size();k++) {
      size_t i = before.XOffset+k;

      Result_Variable_t var;
      const string& varname = APLCON_::BuildVarName(name, before.Values.size(), k);
      var.PristineName = name;
      var.Dimension = before.Values.size();
      var.Index = k;
      var.Value = {*(before.Values[k]), X[i]};

      // sigma is sqrt of diagonal element in V
      const size_t V_ii = APLCON_::V_ij(i,i);
      const double after_sigma = sqrt(V[V_ii]);
      var.Sigma = {*(before.Sigmas[k]), after_sigma};

      // only copy stuff back if variable is not internally stored
      // which is indicated by an empty internal store
      if(before.StoredValues.empty())
        *(before.Values[k]) = X[i];
      if(before.StoredSigmas.empty())
        *(before.Sigmas[k]) = after_sigma;
      if(!before.Pulls.empty())
        *(before.Pulls[k]) = pulls[i];

      // pulls / settings
      var.Pull = pulls[i];
      var.Settings = before.Settings[k];

      // iterating over variables should be the right order
      result.Variables[varname] = var;
    }
  }

  if(fit_settings.SkipCovariancesInResult)
      return result;

  // build the covariances for each variable
  // use the symmetry of V to make it as fast as possible
  for(auto it_zipped_i : APLCON_::index(result.Variables)) {

    for(auto it_zipped_j : APLCON_::index(result.Variables)) {
      const size_t j = it_zipped_j.first;
      const size_t i = it_zipped_i.first;

      if(i>j)
        continue;

      auto& it_map_i = it_zipped_i.second;
      auto& it_map_j = it_zipped_j.second;

      const string& varname_i = it_map_i.first;
      Result_Variable_t& var_i = it_map_i.second;
      const string& varname_j = it_map_j.first;
      Result_Variable_t& var_j = it_map_j.second;

      const size_t V_ij = APLCON_::V_ij(i,j);

      // we use hinted insertion which improves performance by 10%
      // correlations can be obtained

      var_i.Covariances.Before.insert(var_i.Covariances.Before.end(),
                                      make_pair(varname_j, V_before[V_ij]));
      var_i.Covariances.After .insert(var_i.Covariances.After.end(),
                                      make_pair(varname_j, V[V_ij]));
      if(i == j)
        continue;

      // note that V_ij = V_ji

      var_j.Covariances.Before.insert(var_j.Covariances.Before.end(),
                                      make_pair(varname_i, V_before[V_ij]));
      var_j.Covariances.After .insert(var_j.Covariances.After.end(),
                                      make_pair(varname_i, V[V_ij]));
    }
  }



  // consider linked covariances
  for(const auto& it_map : covariances) {
    const covariance_t& cov = it_map.second;
    // don't copy back internally stored values
    if(!cov.StoredValues.empty())
      continue;
    for(size_t i=0;i<cov.Values.size();i++) {
      double* p = cov.Values[i];
      // not all covariances might be linked
      if(p==nullptr)
        continue;
      *p = V[cov.V_ij[i]];
    }
  }

  return result;
}


void APLCON::Init()
{
  // check if we can do some quick init
  if(initialized && instance_id == instance_lastfit) {
    // fully init APLCON,
    InitAPLCON();

    // reset the covariance matrix
    V = V_before;

    // copy again the linked variables to X
    for(const auto& it_map : variables) {
      const variable_t& var = it_map.second;
      auto X_offset = X.begin() + var.XOffset;
      auto dereference = [] (const double* d) {return *d;};
      transform(var.Values.begin(), var.Values.end(), X_offset, dereference);
      // copy the sigmas to diagonal of V (with additional square operation)
      APLCON_::V_transform(V, var.Sigmas, var.V_ij,
                              [] (double d) {return pow(d,2);});
    }

    // copy again the true non-diagonal covariances
    // the distinction between sigmas and covariances makes the interface hopefully more usable,
    // because always specifying covariances is tedious, but sigmas are crucial
    // for measured variables or, say, constrained fitting
    for(const auto& it_map : covariances) {
      const auto& cov = it_map.second;
      APLCON_::V_transform(V, cov.Values, cov.V_ij);
    }
    return;
  }


  // build the storage arrays X, V, F for APLCON

  // X are simply the start values, but also track the
  // map of variables names to index in X (as offsets)
  // this is used to create the double pointer arrays
  // for the constraints later and also to unmap results of APLCON in DoFit
  nVariables = 0;
  X.clear();
  V.clear();
  for(auto& it_map : variables) {
    variable_t& var = it_map.second;
    size_t offset = X.size();
    var.XOffset = offset;

    // in Init is the best place to do this pointer business,
    // so make it for the internally stored variables here
    // we need to work on Values/Sigmas independently, since Sigmas are only optionally linked
    APLCON_::make_pointers_if_any(var.StoredValues, var.Values);
    APLCON_::make_pointers_if_any(var.StoredSigmas, var.Sigmas);

    // now, externally linked variables and internally stored can
    // be treated equally
    nVariables += var.Values.size();
    var.V_ij.resize(var.Values.size());
    for(size_t i=0;i<var.Values.size();i++) {
      X.push_back(*(var.Values[i])); // copy initial values to X

      // take care of diagonal elements in V,
      // off diagonal elements are set to zero by resize (see also covariance handling later)
      // this is like a push_back but with some padding in between
      const size_t j = offset+i;
      const size_t V_ij = APLCON_::V_ij(j,j);
      var.V_ij[i] = V_ij; // remember for later (see covariance init below)
      V.resize(V_ij+1, 0);
      V.back() = pow(*(var.Sigmas[i]),2); // last element is sigma^2
    }
  }



  // F will be set by APLCON iteration loop in DoFit
  // F_func are bound to the double pointers which we know
  // since X is now finally allocated in memory
  nConstraints = 0;
  F_func.clear();
  F_func.reserve(constraints.size());
  for(auto& it_map : constraints) {
    // build the vector of double pointers
    constraint_t& constraint = it_map.second;
    vector< vector<const double*> > args;
    args.reserve(constraint.VariableNames.size()); // args usually smaller, but probably not larger (but not excluded)
    for(const string& varname : constraint.VariableNames) {
      const variable_t& var = GetVariableByName(
            varname,
            "Constraint '"+it_map.first+"' refers to unknown variable '"+varname+"'");
      // check if constraint fits to variables
      if(constraint.WantsDouble && var.Values.size()>1) {
        stringstream msg;
        msg << "Constraint '" << it_map.first << "' wants only single double arguments, "
            << "but '" << varname << "' consists of " << var.Values.size() << " (i.e. more than 1) values.";
        throw Error(msg.str());
      }
      // build the vector of pointers to X values
      vector<const double*> p(var.Values.size());
      const auto& X_offset = X.begin()+var.XOffset;
      transform(X_offset, X_offset+p.size(), p.begin(), APLCON_::make_pointer<double>);
      // and store it in args
      args.push_back(p);
    }
    const auto& func = bind(constraint.Function, args);
    // now, since we have bound the func, we can execute it once
    // to determine the returned number of values and
    // thus obtain the number of constraints
    const vector<double>&  r = func();
    nConstraints += r.size();
    constraint.Number = r.size();
    F_func.push_back(func);
  }
  F.resize(nConstraints);

  // V filled with off-diagonal elements from covariances
  // the variables already have their Values pointer correctly filled,
  // so it's size can be used for the variables's dimensionality
  // V is composed of submatrices due to this different dimensionality.
  for(auto& it_map : covariances) {
    const pair<string, string>& varnames = it_map.first;
    covariance_t& cov = it_map.second;

    // create the pointers for internally stored covariances now,
    // because we need cov.Values.size() to be correct in the following
    // no matter if linked of internal covariance is used
    APLCON_::make_pointers_if_any(cov.StoredValues, cov.Values);

    const string& cov_name = "<'"+varnames.first+"','"+varnames.second+"'>";

    // the setup of the index mapping filed V_ij in constraint_t
    // differs somewhat for off-diagonal vs. diagonal covariance elements
    // so we control this with a little flag
    const bool varnames_equal = varnames.first == varnames.second;

    const variable_t& var1 = GetVariableByName(
          varnames.first,
          "Variable name '"+varnames.first+"' for covariance "+cov_name+" not defined");

    // only search second varname if unequal
    const variable_t& var2 = varnames_equal
        ? var1 :
          GetVariableByName(
                   varnames.second,
                   "Variable name '"+varnames.second+"' for covariance "+cov_name+" not defined");

    const size_t n1 = var1.Values.size();
    const size_t n2 = var2.Values.size();

    if(varnames_equal && n1==1) {
      throw Error("Use sigma to define uncertainty of scalar covariance "+cov_name);
    }

    // expected size of the submatrix without diagonal elements (if varnames equal)
    const size_t v_n = varnames_equal ? n1*(n1-1)/2 : n1*n2;

    if(v_n != cov.Values.size()) {
      stringstream msg;
      msg << "Covariance " << cov_name << " provides " << cov.Values.size()
          << " element" << (cov.Values.size()==1?"":"s") << ", but " << v_n << " covariances needed with"
          << " variable dimensions <" << n1 << "," << n2 << ">";
      throw Error(msg.str());
    }




    // build the indices V_ij, depending on XOffsets of the variables
    // var1 refers to rows, var2 to columns (standard mathematics convention)
    cov.V_ij.reserve(cov.Values.size());
    for(size_t i=0;i<n1;i++) {
      for(size_t j=0;j<n2;j++) {

        // again, handle the special case when varnames are equal
        //const size_t i_ = varnames_equal ?
        if(varnames_equal && i<=j)
          continue;

        // also, check if covariance defined for unmeasured variable
        // which is not meaningful, I guess
        const double s1 = *(var1.Sigmas[i]);
        const double s2 = *(var2.Sigmas[j]);
        // calculating the position is different for diagonal/off-diagonal
        // use V_ij as offset of (i-1) x (i-1) large matrix,
        // note that i=j=0 is excluded due to i<j condition
        const size_t v_ij = varnames_equal ? APLCON_::V_ij(i-1,j) : i*n2+j;
        // make sure cov.Values entry p is valid,
        // then see if this covariance connects unmeasured variables
        const double* p = cov.Values[v_ij];
        if(APLCON_::V_validentry(p) &&
           (s1 == 0 || s2 == 0)
           ) {
          // valid cov entry, but at least one variable is set to "unmeasured"
          // figure out which one to provide helpful error message
          stringstream ss;
          ss << "Variable";
          if(s1 == 0 && s2 != 0) {
            ss << " " << APLCON_::BuildVarName(varnames.first,  n1, i);
          }
          else if(s1 != 0 && s2 == 0) {
            ss << " " << APLCON_::BuildVarName(varnames.second, n2, j);
          }
          else {
            ss << "s "
               << APLCON_::BuildVarName(varnames.first,  n1, j)
               << APLCON_::BuildVarName(varnames.second, n2, j);
          }
          ss << " in covariance "+cov_name+ " has vanishing sigma, i.e. is unmeasured";
          throw Error(ss.str());
        }

        // V_ij with offsets from corresponding variables
        const size_t V_ij = APLCON_::V_ij(var1.XOffset+i,var2.XOffset+j);
        cov.V_ij.push_back(V_ij);
      }
    }

    // now, with some properly initialized cov.V_ij for each case,
    // we can fill the non-diagonal values of V
    APLCON_::V_transform(V, cov.Values, cov.V_ij);
  }

  // finally we know the number of variables and constraints
  // so we can setup APLCON itself
  InitAPLCON();


  // save a pristine copy for later
  V_before = V;

  // remember that this instance has inited APLCON
  initialized = true;
  instance_lastfit = instance_id;
}

void APLCON::InitAPLCON() {

  c_aplcon_aplcon(nVariables, nConstraints);

  c_aplcon_aprint(6, fit_settings.DebugLevel); // default output on LUNP 6 (STDOUT)
  if(isfinite(fit_settings.ConstraintAccuracy))
    c_aplcon_apdeps(fit_settings.ConstraintAccuracy);
  if(isfinite(fit_settings.Chi2Accuracy))
    c_aplcon_apepschi(fit_settings.Chi2Accuracy);
  if(fit_settings.MaxIterations>=0)
    c_aplcon_apiter(fit_settings.MaxIterations);
  if(isfinite(fit_settings.MeasuredStepSizeFactor))
    c_aplcon_apderf(fit_settings.MeasuredStepSizeFactor);
  if(isfinite(fit_settings.UnmeasuredStepSizeFactor))
    c_aplcon_apderu(fit_settings.UnmeasuredStepSizeFactor);
  if(isfinite(fit_settings.MinimalStepSizeFactor))
    c_aplcon_apdlow(fit_settings.MinimalStepSizeFactor);

  for(const auto& it_var : variables) {
    const variable_t& var = it_var.second;
    for(size_t j=0;j<var.Settings.size();j++) {
      const Variable_Settings_t& s = var.Settings[j];
      const int i = j+var.XOffset+1; // APLCON/Fortran starts counting at 1
      // setup APLCON variable specific things
      switch (s.Distribution) {
      case APLCON::Distribution_t::Gaussian:
        // thats the APLCON default, nothing must be called
        break;
      case APLCON::Distribution_t::Poissonian:
        c_aplcon_apoiss(i);
        break;
      case APLCON::Distribution_t::LogNormal:
        c_aplcon_aplogn(i);
        break;
      case APLCON::Distribution_t::SquareRoot:
        c_aplcon_apsqrt(i);
        break;
        // APLCON exposes even more transformations (see wrapper),
        // but they're not mentioned in the README...
      default:
        break;
      }

      if(isfinite(s.Limit.Low) && isfinite(s.Limit.High))
        c_aplcon_aplimt(i, s.Limit.Low, s.Limit.High);
      if(isfinite(s.StepSize))
        c_aplcon_apstep(i, s.StepSize);
    }
  }
}
