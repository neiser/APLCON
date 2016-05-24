#ifndef _APLCON_APLCON_CC_HPP
#define _APLCON_APLCON_CC_HPP 1

#include <algorithm>
#include <functional>
#include <sstream>
#include <vector>
#include <cmath>



namespace APLCON_ {

template<typename T>
T* make_pointer(T& t) {
  return std::addressof(t);
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

// define our own zipped iterator which has an internal index counter
// stolen from http://stackoverflow.com/a/10962575

template <typename T>
struct iterator_extractor { typedef typename T::iterator type; };

template <typename T>
struct iterator_extractor<T const> { typedef typename T::const_iterator type; };


template <typename T>
class Indexer {
public:
    class iterator {
        typedef typename iterator_extractor<T>::type inner_iterator;

        typedef typename std::iterator_traits<inner_iterator>::reference inner_reference;
    public:
        typedef std::pair<size_t, inner_reference> reference;

        iterator(inner_iterator it): _pos(0), _it(it) {}

        reference operator*() const { return reference(_pos, *_it); }

        iterator& operator++() { ++_pos; ++_it; return *this; }
        iterator operator++(int) { iterator tmp(*this); ++*this; return tmp; }

        bool operator==(iterator const& it) const { return _it == it._it; }
        bool operator!=(iterator const& it) const { return !(*this == it); }

    private:
        size_t _pos;
        inner_iterator _it;
    };

    Indexer(T& t): _container(t) {}

    iterator begin() const { return iterator(_container.begin()); }
    iterator end() const { return iterator(_container.end()); }

private:
    T& _container;
}; // class Indexer

template <typename T>
Indexer<T> index(T& t) { return Indexer<T>(t); }


} // end namespace APLCON_

#endif

