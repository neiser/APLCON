#include <iostream>

using namespace std;


#include <APLCON.hpp>

class A {

};

int main() {


  APLCON a("test");


  //cout << a.DoFit() << endl;

  struct BF_t {
    double A;
    double B;
  };

  BF_t BF_e   = {0.105,0.135};
  BF_t BF_tau = {0.095,0.14};

  auto linker = [] (BF_t& v) -> vector<double*> { return {&v.A, &v.B}; };

  a.LinkVariable("BF_e", linker(BF_e), {0.01,0.03});
  a.LinkVariable("BF_tau", linker(BF_tau), {0.03,0.03});

//  auto vector_equal = [] (vector<const double*> a, vector<const double*> b) -> vector<double> {
//    // TODO: check check if sizes of a and b are equal
//    vector<double> r(a.size());
//    for(size_t i=0;i<a.size();i++)
//      r[i] = *a[i]-*b[i];
//    return r;
//  };

  auto make_equal = [] (vector<const double*> a) -> vector<double> {
    return {*a[0] - *a[1]};
  };

  a.AddConstraint("BF_e_equal",{"BF_e"},make_equal);
  a.AddConstraint("BF_tau_equal",{"BF_tau"},make_equal);

  cout << a.DoFit() << endl;


//  struct Vec {
//    double E;
//    double px;
//    double py;
//    double pz;
//  };

//  Vec vec1 = {0,0,0,0};
//  Vec vec2 = {0,0,0,0};
//  auto linker_E = [] (Vec& v) -> vector<double*> { return {&v.E}; };
//  auto linker_p = [] (Vec& v) -> vector<double*> { return {&v.px, &v.py, &v.pz}; };

//  a.LinkVariable("Vec1_E", linker_E(vec1), {1});
//  a.LinkVariable("Vec1_p", linker_p(vec1), {2});
//  a.LinkVariable("Vec2_E", linker_E(vec2), {3});
//  a.LinkVariable("Vec2_p", linker_p(vec2), {3});

//  auto equal_energy = [] (vector<const double*> a, vector<const double*> b) -> vector<double> {
//    return {*a[0] - *b[0]}; // might return more than one constraint
//  };

//  auto equal_momentum_3 = [] (vector<const double*> a, vector<const double*> b) -> vector<double> {
//    // one should check that the vectors a, b have the appropiate lengths...
//    return {
//      *a[0] - *b[0],
//      *a[1] - *b[1],
//      *a[2] - *b[2]
//    }; // might return more than one constraint
//  };

//  auto vector_equal = [] (vector<const double*> a, vector<const double*> b) -> vector<double> {
//    // TODO: check check if sizes of a and b are equal
//    vector<double> r(a.size());
//    for(size_t i=0;i<a.size();i++)
//      r[i] = *a[i]-*b[i];
//    return r;
//  };

//  a.AddConstraint("equal_energy",    {"Vec1_E", "Vec2_E"}, equal_energy);
//  a.AddConstraint("equal_momentum",  {"Vec1_p", "Vec2_p"}, equal_momentum_3);
//  //a.AddLinkedConstraint<3>("equal_momentum_", {"Vec1_E", "Vec2_E"}, equal_momentum_3);

//  a.DoFit();

  //cout << vec1.E << endl;
}
