#include <iostream>
#include <APLCON.hpp>
#include <cmath>
#include <iomanip>

using namespace std;

int main() {

  APLCON a("Example");
  a.AddMeasuredVariable("BF_e_A",   0.1050, 0.01);
  a.AddMeasuredVariable("BF_e_B",   0.135,  0.03);
  a.AddMeasuredVariable("BF_tau_A", 0.095,  0.03);
  a.AddMeasuredVariable("BF_tau_B", 0.14,   0.03);

  auto make_equal = [] (double a, double b) { return a - b; };
  a.AddConstraint("BF_e_equal", {"BF_e_A", "BF_e_B"}, make_equal);
  a.AddConstraint("BF_tau_equal", {"BF_tau_A", "BF_tau_B"}, make_equal);

  //a.AddCovariance("BF_e_A", "BF_tau_B", 0.1);
  //a.AddCovariance("BF_e_B", "BF_tau_A", 0.2);

  const APLCON::Result_t& r = a.DoFit();
  cout << r << endl;

  //APLCON b(a, "Another example");
  //a.AddConstraint("BF_e_equal_", {"BF_e_A", "BF_e_B"}, make_equal);
  //cout << b.DoFit() << endl;
  cout << a.DoFit() << endl;

}
