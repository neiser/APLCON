#include <iostream>
#include <APLCON.hpp>

using namespace std;

//vector<double> make_equal_(const vector<const double*>& a) {
//  return  {a[0] - a[1]};
//}

int main() {

  APLCON a("Linker");

  struct BF_t {
    double A;
    double B;
  };

//  BF_t BF_e   = {0.105,0.135};
//  BF_t BF_tau = {0.095,0.14};

//  auto linker = [] (BF_t& v) -> vector<double*> { return {&v.A, &v.B}; };

  //a.LinkVariable("BF_e", linker(BF_e), {0.01,0.03});
  //a.LinkVariable("BF_tau", linker(BF_tau), {0.03,0.03});

  auto scalar2scalar = [] (const double& a) -> double {
    return  a;
  };

  auto scalar2vector = [] (const double& a) -> vector<double> {
    return  {a};
  };

  auto vector2scalar = [] (const vector<double>& a) -> double {
    return  a[0] - a[1];
  };

  auto vector2vector = [] (const vector<double>& a) -> vector<double> {
    return  {a[0] - a[1]};
  };

  a.AddConstraint("1", {"BF_e"}, scalar2scalar);
  a.AddConstraint("2", {"BF_e"}, scalar2vector);
  a.AddConstraint("3", {"BF_e"}, vector2scalar);
  a.AddConstraint("4", {"BF_e"}, vector2vector);

  auto wrong = [] (const vector<double>& a, double b) -> vector<double> {
    return  {a[0] - a[1] + b};
  };

  //a.AddConstraint("4", {"BF_e", "BF_tau"}, wrong);


  //a.AddConstraint("BF_e_equal",   {"BF_e"}, make_equal);
  //a.AddConstraint("BF_tau_equal", {"BF_tau"}, make_equal);

  //cout << a.DoFit() << endl;

  //cout << "Updated value: " << BF_e.A << endl;
}
