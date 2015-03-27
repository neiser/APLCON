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

  // in this example, the values stay the same,
  // but covariances and the lepton-universality constraint are changed
  // please note that it's not good to switch between a/b APLCON instances
  // during fits, since the APLCON backend needs to initialized every time

  APLCON a("WITHOUT lepton universality");

  // variables are identified by strings,
  // value and sigma specified for the measured variables
  // (you optionally specify distributions, and so on)
  a.AddMeasuredVariable("BF_e_A",   0.1050, 0.01);
  a.AddMeasuredVariable("BF_e_B",   0.135,  0.03);
  a.AddMeasuredVariable("BF_tau_A", 0.095,  0.03);
  a.AddMeasuredVariable("BF_tau_B", 0.14,   0.03);

  // setup a lambda function, which returns 0 if the provided parameters are equal
  auto make_equal = [] (double a, double b) { return a - b; };

  // those two constraints (fulfilled when returning zero)
  // require the equality of the two measured values
  // for the BF of each lepton
  // this is achieved by mapping the lambda function to the named variables
  a.AddConstraint("BF_e_equal", {"BF_e_A", "BF_e_B"}, make_equal);
  a.AddConstraint("BF_tau_equal", {"BF_tau_A", "BF_tau_B"}, make_equal);

  // copy a to b with different name, and add another constraint
  APLCON b(a, "WITH lepton universality");
  b.AddConstraint("BF_equal", {"BF_e_A", "BF_tau_A"}, make_equal);

  // now we can finally do the different fits (some with correlations)

  // first, no correlation
  // note that DoFit() returns a pretty comprehensive struct of type APLCON::Result_t
  // using the variable also makes the APLCON banner appear at the right place, by the way :)
  const APLCON::Result_t& r = a.DoFit();
  cout << "== No Correlation" << endl;
  cout << r << endl;
  cout << b.DoFit() << endl;

  // positive correlation for same observable
  a.SetCovariance("BF_e_A", "BF_e_B", 0.45e-4);
  b.SetCovariance("BF_e_A", "BF_e_B", 0.45e-4);
  // then fit and output
  cout << "== positive correlation, same observable" << endl;
  cout << a.DoFit() << endl;
  cout << b.DoFit() << endl;

  // negative correlation for same observable
  a.SetCovariance("BF_e_A", "BF_e_B", -0.45e-4);
  b.SetCovariance("BF_e_A", "BF_e_B", -0.45e-4);
  // then fit and output
  cout << "== negative correlation, same observable" << endl;
  cout << a.DoFit() << endl;
  cout << b.DoFit() << endl;

  // since we're using the same instances a/b for all four cases,
  // don't forget to clear the covariances for the next examples
  a.SetCovariance("BF_e_A", "BF_e_B", 0);
  b.SetCovariance("BF_e_A", "BF_e_B", 0);

  // positive correlation for different observable
  a.SetCovariance("BF_e_B", "BF_tau_B", 8.96e-4);
  b.SetCovariance("BF_e_B", "BF_tau_B", 8.96e-4);
  // then fit and output
  cout << "== positive correlation, different observable" << endl;
  cout << a.DoFit() << endl;
  cout << b.DoFit() << endl;

  // negative correlation for different observable
  a.SetCovariance("BF_e_B", "BF_tau_B", -8.96e-4);
  b.SetCovariance("BF_e_B", "BF_tau_B", -8.96e-4);
  // then fit and output
  cout << "== negative correlation, different observable" << endl;
  cout << a.DoFit() << endl;
  cout << b.DoFit() << endl;

  return 0;
}
