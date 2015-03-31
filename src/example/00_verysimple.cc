#include <iostream>
#include <APLCON.hpp>

using namespace std;

int main() {

  // this example shows how to do standard Gaussian error propagation
  // with two measured variables A and B
  // APLCON will calculate their sum C=A+B and the propagated error

  APLCON a("Error propagation");

  a.AddMeasuredVariable("A", 10, 0.3);
  a.AddMeasuredVariable("B", 20, 0.4);

  a.AddUnmeasuredVariable("C"); // default value 0, unmeasured means sigma=0

  // setup a lambda function which returns 0
  // if C=A+B aka C - A - B = 0 holds
  auto equality_constraint = [] (double a, double b, double c) { return c - a - b; };
  a.AddConstraint("A+B=C", {"A", "B", "C"}, equality_constraint);

  // do the fit, obtain ra structure
  const APLCON::Result_t& ra = a.DoFit();
  cout << ra << endl;

  // this shows what can access in the result structure ra
  cout << "C's value (should be 30 due to constraint):         "
       << ra.Variables.at("C").Value.After << endl;
  cout << "C's sigma (should be 0.5 due to error propagation): "
       << ra.Variables.at("C").Sigma.After << endl;
  cout << "Correlation between C and B:                        ";
  cout << 100*ra.Variables.at("C").Correlations.After.at("B") << " %" << endl << endl;

  // let's try the same with Poissonian variables
  APLCON b("Poissonian error propagation");

  APLCON::Variable_Settings_t settings = APLCON::Variable_Settings_t::Default;
  settings.Distribution = APLCON::Distribution_t::Poissonian;

  b.AddMeasuredVariable("A", 10, 1, settings);
  b.AddMeasuredVariable("B", 20, 2, settings);

  b.AddUnmeasuredVariable("C");

  b.AddConstraint("A+B=C", {"A", "B", "C"}, equality_constraint);

  const APLCON::Result_t& rb = b.DoFit();
  cout << rb << endl;

  return 0;
}
