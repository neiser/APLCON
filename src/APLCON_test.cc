#include <iostream>

using namespace std;


#include <APLCON.hpp>

class A {

};

int main() {


  APLCON a("test");

  a.AddMeasuredVariable("BF_e_A",   0.1050, 0.01);
  a.AddMeasuredVariable("BF_e_B",   0.135,  0.03);
  a.AddMeasuredVariable("BF_tau_A", 0.095,  0.03);
  a.AddMeasuredVariable("BF_tau_B", 0.14,   0.03);

  auto make_equal = [] (double a, double b) { return a - b; };
  a.AddConstraint("BF_e_equal", {"BF_e_A", "BF_e_B"}, make_equal);
  a.AddConstraint("BF_tau_equal", {"BF_tau_A", "BF_tau_B"}, make_equal);

  //cout << a.DoFit() << endl;


  struct Vec {
    double E;
    double px;
    double py;
    double pz;
  };

  Vec vec1 = {0,0,0,0};
  Vec vec2 = {0,0,0,0};
  auto linker = [] (Vec& v) -> vector<double*> { return {&v.E}; };

  a.LinkVariable("Vec1", vec1, linker);
  a.LinkVariable("Vec2", vec2, linker);

  a.Test("Vec1");
  a.Test("Vec1");
  vec1.E = 5;
  a.Test("Vec1");

  cout << &vec1.E << endl;
  vec1 = vec2;
  cout << &vec1.E << endl;

  a.Test("Vec1");
  vec1.E = 6;
  a.Test("Vec1");
  a.Test("Vec1");


  auto equal_energy = [] (vector<const double*> a, vector<const double*> b) -> vector<double> {
    return {*a[0] - *b[0]}; // might return more than one constraint...?
  };

  a.AddLinkedConstraint("equal_energy", {"Vec1", "Vec2"}, equal_energy);

  //cout << vec1.E << endl;
}
