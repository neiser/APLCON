#ifndef _APLCON_APLCON_OSTREAM_HPP
#define _APLCON_APLCON_OSTREAM_HPP 1

#include "APLCON.hpp"

#include <cmath>
#include <string>
#include <ostream>
#include <iomanip>
#include <iterator>
#include <stdexcept>
#include <sstream>
#include <vector>


// some defaults for PrintFormatting
std::string APLCON::PrintFormatting::Indent = "   ";
std::string APLCON::PrintFormatting::Marker = ">> ";
int APLCON::PrintFormatting::Width = 13;


// some helper string stream which copies the
// formatting "settings" of the given base stream
namespace APLCON_ {

class mystringstream : public std::stringstream {
public:
  mystringstream(const std::ios_base& base) :
    std::stringstream() {
    precision(base.precision());
    flags(base.flags());
  }
};

} // namespace APLCON_

std::ostream& operator<< (std::ostream& o, const APLCON::Limit_t& l) {
  if(!std::isfinite(l.Low) && !std::isfinite(l.High)) {
    return o << "(nolimits)";
  }
  return o << "(" << l.Low << ", " << l.High << ")";
}

std::ostream& operator<< (std::ostream& o, const APLCON::Distribution_t& d) {
  switch (d) {
  case APLCON::Distribution_t::Gaussian:
    o << "Gaussian";
    break;
  case APLCON::Distribution_t::LogNormal:
    o << "LogNormal";
    break;
  case APLCON::Distribution_t::Poissonian:
    o << "Poissonian";
    break;
  case APLCON::Distribution_t::SquareRoot:
    o << "SquareRoot";
    break;
  default:
    throw APLCON::Error("Unkown Distribution_t in ostream");
    break;
  }
  return o;
}

std::ostream& operator<< (std::ostream& o, const APLCON::Result_Status_t& s) {
  switch (s) {
  case APLCON::Result_Status_t::Success:
    o << "Success";
    break;
  case APLCON::Result_Status_t::NoConvergence:
    o << "NoConvergence";
    break;
  case APLCON::Result_Status_t::TooManyIterations:
    o << "TooManyIterations";
    break;
  case APLCON::Result_Status_t::UnphysicalValues:
    o << "UnphysicalValues";
    break;
  case APLCON::Result_Status_t::NegativeDoF:
    o << "NegativeDoF";
    break;
  case APLCON::Result_Status_t::OutOfMemory:
    o << "OutOfMemory";
    break;
  default:
    throw APLCON::Error("Unknown Result_Status_t in ostream");
    break;
  }

  return o;
}

std::ostream& operator<< (std::ostream& o, const APLCON::Variable_Settings_t& s) {
  APLCON_::mystringstream limit(o); // inherit the format properties of o
  limit << s.Limit;
  APLCON_::mystringstream stepsize(o); // inherit the format properties of o
  if(std::isfinite(s.StepSize))
    stepsize << s.StepSize;
  else
    stepsize << "def_stepsize";
  return o << std::left
           << std::setw(11) << s.Distribution // longest distribution name is 10
           << std::setw(15) << limit.str()
           << std::setw(10) << (s.StepSize==0 ? "fixed" : stepsize.str())
           << std::right;
}

template<typename F>
void stringify_covariances(
    std::ostream& o,
    const std::map<std::string, APLCON::Result_Variable_t>& variables,
    const std::string& in,
    F f,
    int w_varname,
    bool skipUnmeasured,
    double factor = 1.0
    ) {

  w_varname += 4;
  const int w = APLCON::PrintFormatting::Width;

  // print top line
  o << in << std::setw(w_varname) << " ";
  for(const auto& it_zipped : APLCON_::index(variables)) {
    const size_t i = it_zipped.first;
    const auto& it_map = it_zipped.second;
    const APLCON::Result_Variable_t& v = it_map.second;
    if(skipUnmeasured && v.Sigma.Before == 0)
      continue;
    std::stringstream i_;
    i_ << "(" << i << ")";
    o << std::setw(w) << i_.str();
  }
  o << std::endl;

  // print matrix

  std::vector<APLCON::Result_Variable_t> v_variables;

  for(const auto& it_zipped : APLCON_::index(variables)) {
    const size_t i = it_zipped.first;
    const auto& it_map = it_zipped.second;
    const std::string varname = it_map.first;
    const APLCON::Result_Variable_t& v = it_map.second;

    v_variables.emplace_back(v);

    if(skipUnmeasured && v.Sigma.Before == 0)
      continue;

    std::stringstream i_;
    const size_t padding = i<10 ? 0 : 1; // breaks with more than 100 variables...
    i_ << i << ") ";
    o << in << std::setw(4-padding) << i_.str()
      << std::left << std::setw(w_varname-4) << varname << std::right;

    const auto& cov = f(varname, v);
    for(const auto& it_zipped : APLCON_::index(cov)) {
      const size_t j = it_zipped.first;
      const auto& it_map = it_zipped.second;
      const double& cov_val = it_map.second;

      if(j>i)
        continue;

      const APLCON::Result_Variable_t& v_j = v_variables[j];
      if(skipUnmeasured && v_j.Sigma.Before == 0)
        continue;
      o << std::setw(w) << cov_val*factor;
    }
    o << std::endl;
  }

  // print bottom line
  o << in << std::setw(w_varname) << " ";
  for(const auto& it_zipped : APLCON_::index(variables)) {
    const size_t i = it_zipped.first;
    const auto& it_map = it_zipped.second;
    const APLCON::Result_Variable_t& v = it_map.second;
    if(skipUnmeasured && v.Sigma.Before == 0)
      continue;
    std::stringstream i_;
    i_ << "(" << i << ")";
    o << std::setw(w) << i_.str();
  }
  o << std::endl;
}

void stringify_variables(
    std::ostream& o,
    const std::map<std::string, APLCON::Result_Variable_t>& variables,
    const std::string& extra_indent = "",
    const bool success = true) {
  // do some extra work and find out the maximum length of
  // variable names
  const int w = APLCON::PrintFormatting::Width;
  size_t w_varname = w;
  for(const auto& it_map : variables) {
    if(w_varname<it_map.first.size())
      w_varname = it_map.first.size();
  }
  w_varname += 2;


  // calculate the correlations
  auto correlations = APLCON::CalculateCorrelations(variables);

  const std::string& in = extra_indent + APLCON::PrintFormatting::Indent;
  const std::string& ma = extra_indent + APLCON::PrintFormatting::Marker;

  // print stuff before the Fit
  o << ma << "Before Fit:" << std::endl << std::endl;
  o << in
    << std::left << std::setw(w_varname) << "Variable" << std::right
    << std::setw(w)   << "Value"
    << std::setw(w)   << "Sigma"
    << std::left      << "   Settings" << std::right
    << std::endl;
  for(const auto& it_map : variables) {
    const APLCON::Result_Variable_t& v = it_map.second;
    APLCON_::mystringstream sigma(o);
    if(v.Sigma.Before != 0)
      sigma << v.Sigma.Before;
    else
      sigma << "unmeas";
    o << in
      << std::left << std::setw(w_varname) << it_map.first << std::right
      << std::setw(w) << v.Value.Before
      << std::setw(w) << sigma.str()
      << std::left << "   " << v.Settings << std::right
      << std::endl;
  }
  o << std::endl;

  if(!correlations.empty()) {
      o << in << "Covariances: " << std::endl;
      stringify_covariances(o, variables, in,
                            [](const std::string&, const APLCON::Result_Variable_t& v) {
          return v.Covariances.Before;
      },
      w_varname,
      true);
      o << std::endl;
      o << in << "Correlations (in %): " << std::endl;
      auto correlations_before = [&correlations]
              (const std::string& varname,
              const APLCON::Result_Variable_t&) {
          return correlations[varname].Before;
      };
      stringify_covariances(o, variables, in,
                            correlations_before,
                            w_varname,
                            true,
                            100);
  }

  if(!success)
    return;

  // print stuff after the fit
  o << ma << "After Fit:" << std::endl << std::endl;
  o << in
    << std::left << std::setw(w) << "Variable" << std::right
    << std::setw(w) << "Value"
    << std::setw(w) << "Sigma"
    << std::setw(w) << "Pull"
    << std::endl;
  for(const auto& it_map : variables) {
    const APLCON::Result_Variable_t& v = it_map.second;
    o << in
      << std::left << std::setw(w) << it_map.first << std::right
      << std::setw(w) << v.Value.After
      << std::setw(w) << v.Sigma.After
      << std::setw(w) << v.Pull
      << std::endl;
  }
  o << std::endl;

  if(!correlations.empty()) {
      o << in << "Covariances: " << std::endl;
      stringify_covariances(o, variables, in,
                            [](const std::string&, const APLCON::Result_Variable_t& v) {
          return v.Covariances.After;
      },
      w_varname,
      false);
      o << std::endl;
      o << in << "Correlations (in %): " << std::endl;
      auto correlations_after = [&correlations]
              (const std::string& varname,
              const APLCON::Result_Variable_t&) {
          return correlations[varname].After;
      };
      stringify_covariances(o, variables, in,
                            correlations_after,
                            w_varname,
                            false,
                            100);

  }

}

std::ostream& operator<< (std::ostream& o, const APLCON::Result_t& r) {
  const std::string& in = APLCON::PrintFormatting::Indent;
  const std::string& ma = APLCON::PrintFormatting::Marker;

  const bool success = r.Status == APLCON::Result_Status_t::Success;

  const std::string& tag = success ? "" : "ERROR ";

  // general info
  o << ma << (r.Name==""?"APLCON":r.Name) << " with " << r.Variables.size() << " variables and "
    << r.NScalarConstraints << " constraints:" << std::endl;
  o << in << tag << r.Status << " after " << r.NIterations << " iterations, " << r.NFunctionCalls << " function calls " << std::endl;
  o << in << "Chi^2 / DoF = " << r.ChiSquare << " / " << r.NDoF << " = " << r.ChiSquare/r.NDoF << std::endl;
  o << in << "Probability = " << r.Probability << std::endl;
  o << in << "Constraints: ";

  for(const auto& it_zipped : APLCON_::index(r.Constraints)) {
    const size_t i = it_zipped.first;
    const auto& it_map = it_zipped.second;
    const APLCON::Result_Constraint_t& c = it_map.second;
    o << it_map.first;
    if(c.Dimension>1)
      o << "[" << c.Dimension << "]";
    if(i < r.Constraints.size()-1)
      o << ", ";
  }
  o << std::endl;

  stringify_variables(o, r.Variables, in, success);
  return o;
}

#endif

