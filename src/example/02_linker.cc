#include <iostream>
#include <iomanip>
#include <APLCON.hpp>

using namespace std;

int main() {
  // we give the first fit from 01_simple.cc as
  // a linked variable example

  // setup your own data structure
  // might be as complicated as you wish,
  // as long as you can get double pointers out of it
  struct BF_t {
    double A;
    double B;
  };

  BF_t BF_e   = {0.105,0.135};
  BF_t BF_tau = {0.095,0.14};

  // now, the fancy part:

  APLCON a("Linked variables");

  // ah snap, that wasn't fancy, now here we are:

  // we use LinkVariable to give vector<double*> to our own data structures
  // of course, YOU MUST make sure that your data structures are not destroyed
  // when calling the DoFit function

  // so we define some handy little function which converts instances of BF_t
  // note that linker needs an explicit return type declaration " -> vector<double*> "
  // since we use an initializer_list { }
  auto linker = [] (BF_t& v) -> vector<double*> { return {&v.A, &v.B}; };

  a.LinkVariable("BF_tau",       // name
                 linker(BF_tau), // give double* vector
                 // give sigmas for field A and B,
                 // {0.03} would also work (gets expanded for all values)
                 // on newer gcc versions, you may even
                 // drop the vector<double> declaraton...
                 vector<double>{0.03,0.03}
                 );
  // second variable. You can optionally specify settings,
  // and if you provide a sigma of 0 then the
  a.LinkVariable("BF_e", linker(BF_e), vector<double>{0.01,0.03});

  // accordingly, when defining constraints, you get vectors as arguments
  // in the order the linker function (from above) returns it
  auto equality_constraint = [] (vector<double> a) { return a[0] - a[1]; };
  a.AddConstraint("BF_e_equal",   {"BF_e"},   equality_constraint);
  a.AddConstraint("BF_tau_equal", {"BF_tau"}, equality_constraint);

  // do fit and print everything
  // that is identical to the first example,
  // even though we specified the variables in a different order
  cout.precision(3); // set precision globally, which makes output nicer
  cout << a.DoFit() << endl;

  // double-check if linking actually works
  cout << "+++++ Value of BF_e[0]: " << BF_e.A << endl << endl;


  // but there's more when linking variables. we can also link sigmas
  // (full covariances not implemented yet)

  APLCON b("Linked sigmas");

  // we need some data structure to hold it
  BF_e   = {0.105, 0.135};
  BF_tau = {0.095, 0.14};
  BF_t BF_e_sigma = {0.01, 0.03};

  // this is an example, so let's make it easy and reuse the little linker
  // and the constraint from above
  // Note: reusing/modifying the instance "a" is considered stupid, that's why there's
  // no method to re-link variables. One fitter instance should serve one purpose.

  // second variable. You can optionally specify settings,
  // and if you provide a sigma of 0 then the
  b.LinkVariable("BF_e",   linker(BF_e),   linker(BF_e_sigma));
  b.LinkVariable("BF_tau", linker(BF_tau), vector<double>{0.03,0.03});

  b.AddConstraint("BF_e_equal",   {"BF_e"},   equality_constraint);
  b.AddConstraint("BF_tau_equal", {"BF_tau"}, equality_constraint);

  cout << b.DoFit() << endl;

  // double-check if linking actually works
  cout << "+++++ Value of BF_e[0]: " << BF_e.A       << endl;
  cout << "+++++ Sigma of BF_e[0]: " << BF_e_sigma.A << endl << endl;

  // feel free to play around, hopefully you get error messages when you need them,
  // and correct results when you expect them :)
}
