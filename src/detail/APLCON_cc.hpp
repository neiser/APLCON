#ifndef _APLCON_APLCON_CC_HPP
#define _APLCON_APLCON_CC_HPP 1

#include <type_traits>
#include <vector>
#include <functional>
#include <algorithm>

namespace APLCON_ {

template<typename T>
T* make_pointer(T& t) {
  return &t;
}

template <typename T>
void copy_pointers(std::vector<T>& source, std::vector<T*>& target)
{
  if(source.size()==0)
    return;
  target.resize(source.size());
  std::transform(source.begin(), source.end(), target.begin(), make_pointer<T>);
}

} // end namespace APLCON_

#endif

