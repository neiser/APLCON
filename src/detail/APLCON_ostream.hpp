#ifndef _APLCON_APLCON_OSTREAM_HPP
#define _APLCON_APLCON_OSTREAM_HPP 1

#include "APLCON.hpp"

#include <cmath>
#include <string>
#include <ostream>
#include <iomanip>
#include <stdexcept>
#include <sstream>
#include <vector>


// some defaults for PrintFormatting
std::string APLCON::PrintFormatting::Indent = "   ";
std::string APLCON::PrintFormatting::Marker = ">> ";
int APLCON::PrintFormatting::Width = 13;


// some helper string stream which copies the 
// formatting "settings" of the given base stream
class mystringstream : public std::stringstream {
public:
  mystringstream(const std::ios_base& base) :
    std::stringstream() {
    precision(base.precision());
    flags(base.flags());
  }
};

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
  mystringstream limit(o); // inherit the format properties of o 
  limit << s.Limit;
  mystringstream stepsize(o); // inherit the format properties of o 
  if(std::isfinite(s.StepSize))
    stepsize << s.StepSize;
  return o << std::left 
           << std::setw(11) << s.Distribution // longest distribution name is 10
           << std::setw(15) << limit.str()    
           << std::setw(10) << (s.StepSize==0 ? "fixed" : stepsize.str()) 
           << std::right;
}

template<typename F>
void stringify_covariances(
    std::ostream& o,
    const std::vector<APLCON::Result_Variable_t>& variables,
    const std::string& in,
    F f,
    bool skipUnmeasured,
    double factor = 1.0
    ) {

  const int w = APLCON::PrintFormatting::Width;
  const int w_varname = APLCON::PrintFormatting::Width;
  o << in << std::setw(w_varname) << " ";
  for(size_t i=0;i<variables.size();i++) {
    const APLCON::Result_Variable_t& v = variables[i];
    if(skipUnmeasured && v.Sigma.Before == 0)
      continue;
    std::stringstream i_;
    i_ << "(" << i << ")";
    o << std::setw(w) << i_.str();
  }
  o << std::endl;
  for(size_t i=0;i<variables.size();i++) {
    const APLCON::Result_Variable_t& v = variables[i];
    if(skipUnmeasured && v.Sigma.Before == 0)
      continue;
    std::stringstream i_;
    const size_t padding = i<10 ? 0 : 1; // breaks with more than 100 variables...
    i_ << i << ") ";
    o << in << std::setw(4-padding) << i_.str()
      << std::left << std::setw(w_varname-4) << v.Name << std::right;
    const std::vector<double>& cov = f(v);
    for(size_t j=0;j<cov.size();j++) {
      const APLCON::Result_Variable_t& v_j = variables[j];
      if(skipUnmeasured && v_j.Sigma.Before == 0)
        continue;
      o << std::setw(w) << cov[j]*factor;
    }
    o << std::endl;
  }
  o << in << std::setw(w_varname) << " ";
  for(size_t i=0;i<variables.size();i++) {
    const APLCON::Result_Variable_t& v = variables[i];
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
    const std::vector<APLCON::Result_Variable_t>& variables,
    const std::string& extra_indent = "",
    const bool success = true) {
  const int w = APLCON::PrintFormatting::Width;
  const std::string& in = extra_indent + APLCON::PrintFormatting::Indent;
  const std::string& ma = extra_indent + APLCON::PrintFormatting::Marker;
  
  // print stuff before the Fit
  o << ma << "Before Fit:" << std::endl << std::endl;
  o << in
    << std::left << std::setw(w) << "Name" << std::right
    << std::setw(w)   << "Value"
    << std::setw(w)   << "Sigma"
    << std::left      << "   Settings" << std::right
    << std::endl;
  for(const APLCON::Result_Variable_t& v : variables) {
    mystringstream sigma(o);
    if(v.Sigma.Before != 0)
      sigma << v.Sigma.Before;
    else
      sigma << "unmeas";
    o << in
      << std::left << std::setw(w) << v.Name << std::right
      << std::setw(w) << v.Value.Before
      << std::setw(w) << sigma.str()
      << std::left << "   " << v.Settings << std::right
      << std::endl;
  }
  o << std::endl;

  o << in << "Covariances: " << std::endl;
  stringify_covariances(o, variables, in,
                        [](const APLCON::Result_Variable_t& v) {return v.Covariances.Before;},
                        true);

  o << in << "Correlations (in %): " << std::endl;
  stringify_covariances(o, variables, in,
                        [](const APLCON::Result_Variable_t& v) {return v.Correlations.Before;},
                        true,
                        100);

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
  for(const APLCON::Result_Variable_t& v : variables) {
    o << in
      << std::left << std::setw(w) << v.Name << std::right
      << std::setw(w) << v.Value.After
      << std::setw(w) << v.Sigma.After
      << std::setw(w) << v.Pull
      << std::endl;
  }
  o << std::endl;

  o << in << "Covariances: " << std::endl;
  stringify_covariances(o, variables, in,
                        [](const APLCON::Result_Variable_t& v) {return v.Covariances.After;},
                        false);

  o << in << "Correlations (in %): " << std::endl;
  stringify_covariances(o, variables, in,
                        [](const APLCON::Result_Variable_t& v) {return v.Correlations.After;},
                        false,
                        100);

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
  for(size_t i=0;i<r.Constraints.size();i++) {
    const APLCON::Result_Constraint_t& c = r.Constraints[i];
    o << c.Name;
    if(c.Number>1)
      o << "[" << c.Number << "]";
    if(i < r.Constraints.size()-1)
      o << ", ";
  }
  o << std::endl;

  stringify_variables(o, r.Variables, in, success);
  return o;
}

#endif

