#include <iostream>
#include <APLCON.hpp>

using namespace std;

int main() {
  // This is the C++11 version of test/avlass.F in APLCON

  cout << "Combining correlated measurements of several" << endl
       << "different physical quantities" << endl << endl
       << "Fictitious example (A. Valassi, pages 399 - 403)" << endl
       << "A. Valassi, NIMA 500 (203) 391-405" << endl
       << "Case of two experiments A and B, measuring the" << endl
       << "branching fraction of the W boson in the two decay" << endl
       << "channels to electrons and taus." << endl;

  //

  APLCON a("Example");
  a.AddMeasuredVariable("BF_e_A",   0.1050, 0.01);
  a.AddMeasuredVariable("BF_e_B",   0.135,  0.03);
  a.AddMeasuredVariable("BF_tau_A", 0.095,  0.03);
  a.AddMeasuredVariable("BF_tau_B", 0.14,   0.03);

  auto make_equal = [] (double a, double b) { return a - b; };
  a.AddConstraint("BF_e_equal", {"BF_e_A", "BF_e_B"}, make_equal);
  a.AddConstraint("BF_tau_equal", {"BF_tau_A", "BF_tau_B"}, make_equal);

  //a.AddCovariance("BF_e_A", "BF_tau_B", 0.1);

  cout << a.DoFit() << endl;

  //APLCON b(a, "Another example");
  //a.AddConstraint("BF_e_equal_", {"BF_e_A", "BF_e_B"}, make_equal);
  //cout << b.DoFit() << endl;
  //cout << a.DoFit() << endl;

}
