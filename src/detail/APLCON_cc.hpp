#ifndef _APLCON_APLCON_CC_HPP
#define _APLCON_APLCON_CC_HPP 1

#include <algorithm>
#include <functional>
#include <sstream>
#include <vector>



namespace APLCON_ {

template<typename T>
T* make_pointer(T& t) {
  return &t;
}

template <typename T>
void make_pointers_if_any(std::vector<T>& source, std::vector<T*>& target)
{
  if(source.empty())
    return;
  target.resize(source.size());
  std::transform(source.begin(), source.end(), target.begin(), make_pointer<T>);
}

/**
 * @brief V_ij index for symmetric covariance matrix V
 * @param i row index
 * @param j column index
 * @return index for symmetric covariance matrix V
 */
size_t V_ij(const size_t i, const size_t j) {
  // V_ij is symmetric, and APLCON chooses row index <= column index (below diagonal)
  // see also IJSYM in Fortran code, which is identical except that indices start from 0 here
  if(i>j)
    return i*(i+1)/2 + j;
  else
    return j*(j+1)/2 + i;
}

bool V_validentry(const double* p) {
  return p!=nullptr && std::isfinite(*p);
}

void V_transform(
    std::vector<double>& V,
    const std::vector<double*>& values,
    const std::vector<size_t>& V_ij,
    std::function<double(double)> transform = [] (double d) {return d;} // identity function as default
    )
{
  for(size_t i=0;i<values.size();i++) {
    const double* p = values[i];
    // just continue if some value should be kept at default
    if(!V_validentry(p))
      continue;
    V[V_ij[i]] = transform(*p);
  }
}

std::string BuildVarName(const std::string& name, size_t n, size_t k) {
  std::stringstream s_name;
  s_name << name;
  if(n>1) {
    s_name << "[" << k << "]";
  }
  return s_name.str();
}

} // end namespace APLCON_

#endif

