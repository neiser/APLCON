#include <iostream>
#include <APLCON.hpp>

#include "catch.hpp"

using namespace std;

TEST_CASE("Very Simple", "") {

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

  // this shows what can access in the result structure "ra"
  // note that the correlations must be calculated on demand
  cout << "C's value (should be 30 due to constraint):         "
       << ra.Variables.at("C").Value.After << endl;

  REQUIRE(ra.Variables.at("C").Value.After == Approx(30.0));

  cout << "C's sigma (should be 0.5 due to error propagation): "
       << ra.Variables.at("C").Sigma.After << endl;

  REQUIRE(ra.Variables.at("C").Sigma.After == Approx(0.5));

  const auto& correlations = APLCON::CalculateCorrelations(ra.Variables);
  cout << "Correlation between C and B:                        ";
  cout << 100*correlations.at("C").After.at("B") << " %" << endl << endl;

  REQUIRE( correlations.at("C").After.at("B") == Approx(0.8));

}
