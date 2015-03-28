#ifndef _APLCON_APLCON_OSTREAM_HPP
#define _APLCON_APLCON_OSTREAM_HPP 1

#include "APLCON.hpp"

#include <string>
#include <ostream>
#include <iomanip>
#include <stdexcept>
#include <sstream>
#include <vector>


const std::string APLCON::PrintFormatting::Indent = "   ";
const std::string APLCON::PrintFormatting::Marker = ">> ";
const int APLCON::PrintFormatting::Width = 13;

std::ostream& operator<< (std::ostream& o, const APLCON::Limit_t& l) {
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
    throw std::logic_error("Unkown Distribution_t in ostream");
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
    throw std::logic_error("Unkown Result_Status_t in ostream");
    break;
  }

  return o;
}

template<typename F>
std::string stringify_contraints(const std::vector<APLCON::Result_Variable_t>& variables,
                                 const std::string& in,
                                 F f) {
  const int w = APLCON::PrintFormatting::Width;
  const int w_varname = APLCON::PrintFormatting::Width+5;
  std::stringstream o;
  o << in << "Covariances: " << std::endl;
  o << in << std::setw(w_varname) << " ";
  for(size_t i=0;i<variables.size();i++) {
    std::stringstream i_;
    i_ << "(" << i << ")";
    o << std::setw(w) << i_.str();
  }
  o << std::endl;
  for(size_t i=0;i<variables.size();i++) {
    const APLCON::Result_Variable_t& v = variables[i];
    std::stringstream i_;
    i_ << "(" << i << ") ";
    o << in << std::setw(4) << i_.str() << " " << std::setw(w_varname-5) << v.Name;
    const std::vector<double>& cov = f(v);
    for(size_t j=0;j<cov.size();j++) {
      o << std::setw(w) << cov[j];
    }
    o << std::endl;
  }
  o << in << std::setw(w_varname) << " ";
  for(size_t i=0;i<variables.size();i++) {
    std::stringstream i_;
    i_ << "(" << i << ")";
    o << std::setw(w) << i_.str();
  }
  o << std::endl;
  return o.str();
}

std::string stringify_variables(const std::vector<APLCON::Result_Variable_t>& variables,
                                const std::string& extra_indent = "") {
  std::stringstream o;
  const int w = APLCON::PrintFormatting::Width;
  const std::string& in = extra_indent + APLCON::PrintFormatting::Indent;
  const std::string& ma = extra_indent + APLCON::PrintFormatting::Marker;
  // print stuff before the Fit
  o << ma << "Before Fit:" << std::endl << std::endl;
  o << in
    << std::left << std::setw(w) << "Name" << std::right
    << std::setw(w)   << "Value"
    << std::setw(w)   << "Sigma"
    << std::left << std::setw(2*w) << "   Settings" << std::right
    << std::endl;
  for(const APLCON::Result_Variable_t& v : variables) {
    std::stringstream settings;
    settings << "   " << v.Settings.Distribution << " " << v.Settings.Limit << " " << v.Settings.StepSize;
    o << in
      << std::left << std::setw(w) << v.Name << std::right
      << std::setw(w) << v.Value.Before
      << std::setw(w) << v.Sigma.Before
      << std::left << std::setw(2*w) << settings.str() << std::right
      << std::endl;
  }
  o << std::endl;

  o << stringify_contraints(variables, in, [](const APLCON::Result_Variable_t& v) {return v.Covariances.Before;});

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

  o << stringify_contraints(variables, in, [](const APLCON::Result_Variable_t& v) {return v.Covariances.After;});

  return o.str();

}

std::ostream& operator<< (std::ostream& o, const APLCON::Result_t& r) {
  const std::string& in = APLCON::PrintFormatting::Indent;
  const std::string& ma = APLCON::PrintFormatting::Marker;

  // general info
  o << ma << (r.Name==""?"APLCON":r.Name) << " with " << r.Variables.size() << " variables and "
    << r.Constraints.size() << " constraints:" << std::endl;
  o << in << r.Status << " after " << r.NIterations << " iterations, " << r.NFunctionCalls << " function calls " << std::endl;
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

  if(r.Status == APLCON::Result_Status_t::Success)
    o << stringify_variables(r.Variables, in);
  return o;
}

#endif

