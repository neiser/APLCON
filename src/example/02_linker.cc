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


  BF_t BF_e   = {0.105,0.135};
  BF_t BF_tau = {0.095,0.14};

  auto linker = [] (BF_t& v) -> vector<double*> { return {&v.A, &v.B}; };

  a.LinkVariable("BF_e", linker(BF_e), {0.01,0.03});
  a.LinkVariable("BF_tau", linker(BF_tau), {0.03,0.03});

  // this time, we want to have
  auto equality_constraint = [] (double a) { return a - 7; };

  a.AddConstraint("BF_e_equal",   {"BF_e"},   equality_constraint);
  a.AddConstraint("BF_tau_equal", {"BF_tau"}, equality_constraint);

  cout << a.DoFit() << endl;

  //cout << "Updated value: " << BF_e.A << endl;
}
