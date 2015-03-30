#include <iostream>
#include <APLCON.hpp>
#include <limits>
#include <iomanip>

using namespace std;

int main() {

  // this example illustrates some more advanced usage of the APLCON interface,
  // like setting up (optionally linked) covariances
  // and more complex constraint functions
  // it is also a toy model kinematical fitter...

  // please feel free to modify the code here and see
  // if APLCON will throw you a meaningful exception in
  // case you set up something stupid :)


  // this example uses rather stupid values for covariances,
  // so APLCON needs more iterations to converge
  // have a look at APLCON::Fit_Settings_t to see what can be configured globally

  APLCON::Fit_Settings_t settings = APLCON::Fit_Settings_t::Default;
  settings.MaxIterations = 500;
  APLCON a("Fit A", settings);

  // as an example structure, we use something
  // which looks like a Lorentz vector
  struct Vec {
    double E;
    double px;
    double py;
    double pz;
  };

  // just some particles
  // zeroth component should be large enough
  // for non-imaginary invariant mass
  Vec vec1a = { sqrt(4+9+16)*1.02,   2,  3,  4}; // that's a photon up to 2%
  Vec vec2a = { sqrt(4+9+16)*1.05,  -2, -3, -4}; // that's a photon up to 5%, flying in the opposite direction
  Vec vec3a = { 13,                  0,  0,  0}; // that's something with rest mass 13


  // linked sigmas are already demonstrated in 02_linker.cc
  // so here we just assume some errors on the percent level for the photons
  const vector<double> sigma1 = {0.6};
  const vector<double> sigma2 = {0.8};
  // third particle is unmeasured, so sigma=0
  const vector<double> sigma3 = {0};

  // you're totally free in
  // how the fields of your data structure are linked
  // you might have also specified the particle by energy, theta and phi

  // for instance a, we separate E and p
  auto linker_E = [] (Vec& v) -> vector<double*> { return {&v.E}; };
  auto linker_p = [] (Vec& v) -> vector<double*> { return {&v.px, &v.py, &v.pz}; };
  APLCON::Variable_Settings_t fixvar = APLCON::Variable_Settings_t::Default;
  fixvar.StepSize = APLCON::NaN; // set to zero to fix, NaN uses APLCON's default
  a.LinkVariable("Vec1_E", linker_E(vec1a), sigma1);
  a.LinkVariable("Vec1_p", linker_p(vec1a), sigma1);
  a.LinkVariable("Vec2_E", linker_E(vec2a), sigma2, {fixvar});
  a.LinkVariable("Vec2_p", linker_p(vec2a), sigma2);
  a.LinkVariable("Vec3_E", linker_E(vec3a), sigma3);
  a.LinkVariable("Vec3_p", linker_p(vec3a), sigma3);



  // #######################################
  // SetCovariance was already used for scalar variables in the first example
  // here we show how to setup covariances for vector-valued variables


  // EXAMPLE (1): covariances between scalar- and vector-valued variable
  // the first  variable defines the number of rows
  // the second variable defines the number of columns
  // the convention for symmetric covariance matrix is: rows<columns,
  // and variable names according to "first row index, then column index")

  constexpr double E_px = 0.001;
  constexpr double E_py = 0.002;
  constexpr double E_pz = 0.003;
  a.SetCovariance("Vec1_E", "Vec2_p", // variables give 1 row and 3 columns
                  vector<double>{
                    E_px, E_py, E_pz
                  });
  // in this case, swapping the variable names would setup the same covariances,
  // which is not true in general. so always check if you actually made it right :)


  // EXAMPLE (2): covariances for vector-valued variable Vec1_p,
  // which---interpreted as momentum---can have covariances pypx, pzpx, pzpy

  constexpr double pypx = 0.004;
  constexpr double pzpx = 0.005;
  constexpr double pzpy = 0.006;
  // the corresponding covariance matrix is a simple 1-dimensional vector,
  // but interpreted as the part below the diagonal
  // (since the diagonal of the covariance matrix corresponds to sigmas of px, py and pz, which are specified above)
  const vector<double> covariance_matrix_p = {
    // the /**/ indicate the position of the diagonal elements
    /**/             // first  row -> x component
    pypx, /**/       // second row -> y component
    pzpx, pzpy /**/  // third  row -> z component
  };
  a.SetCovariance("Vec1_p", "Vec1_p", covariance_matrix_p);


  // EXAMPLE (3): covariances for two vector-valued variables Vec1_p and Vec2_p
  // since there are no diagonal elements of the full covariance matrix involved,
  // the given vector of covariances represents a 3x3 matrix in this case

  // variables indicate 3 rows and 3 columns
  const vector<double> covariance_matrix_pp = {
   0.001, 0.002, 0.003,
   0.004, 0.005, 0.006,
   0.007, 0.008, 0.009
  };
  a.SetCovariance("Vec1_p", "Vec2_p",covariance_matrix_pp);

  // #######################################
  // now some more advanced constraint business
  // remember that the arguments correspond to the specified
  // in general, there are four different constraint types
  // which are supported by the interface:
  // (1) arguments are all scalars and returns a scalar
  // (2) arguments are all vectors and returns a scalar
  // (3) arguments are all scalars and returns a vector
  // (4) arguments are all vectors and returns a vector
  // in case of (3), (4) the provided constraint aggregates several scalar constraints into one function

  // example for case (2)
  constexpr auto invariant_mass = [] (const vector<double>& E, const vector<double>& p) -> double {
    // note that, although E is a scalar variable,
    // it is provided as a vector with one element
    // (mixing scalar/vector arguments are not supported at the moment)
    // M^2 = E^2 - vec(p)^2
    const double M2 = pow(E[0],2) - pow(p[0],2) - pow(p[1],2) - pow(p[2],2);
    return M2; // require the invariant mass to be zero
  };
  a.AddConstraint("invariant_mass1", {"Vec1_E", "Vec1_p"}, invariant_mass);
  a.AddConstraint("invariant_mass2", {"Vec2_E", "Vec2_p"}, invariant_mass);

  // example for case (4)
  constexpr auto opposite_momentum_3 = [] (const vector<double>& a, const vector<double>& b) -> vector<double> {
    // one may check that the vectors a, b have the appropiate lengths
    // that's something the interface can't do for you...
    return {
      a[0] + b[0],
      a[1] + b[1],
      a[2] + b[2]
    }; // returns 3 scalar constraints
  };
  // the two photons shall fly back to back
  a.AddConstraint("opposite_momentum",  {"Vec1_p", "Vec2_p"}, opposite_momentum_3);

  // to make the fit at least somewhat meaningful, provide the four-momentum conservation,
  // so Vec1+Vec2=Vec3 aka Vec1+Vec2-Vec3 = 0
  constexpr auto require_conservation = [] (
      const vector<double>& v1_E,
      const vector<double>& v1_p,
      const vector<double>& v2_E,
      const vector<double>& v2_p,
      const vector<double>& v3_E,
      const vector<double>& v3_p
      ) -> vector<double> {
    // this is rather tedious to formulate,
    // see instance b below for more elegant solution
    return {
      v1_E[0] + v2_E[0] - v3_E[0],
      v1_p[0] + v2_p[0] - v3_p[0],
      v1_p[1] + v2_p[1] - v3_p[1],
      v1_p[2] + v2_p[2] - v3_p[2]
    };
  };
  a.AddConstraint("require_conservation",
                  {"Vec1_E","Vec1_p","Vec2_E","Vec2_p","Vec3_E","Vec3_p"},
                  require_conservation);

  // don't execute DoFit yet because we copy the initial values of the vectors

  // #######################################
  // we setup instance b exactly as instance a,
  // but this time with 4-vectors
  // this hopefully shows how powerful the interface is :)

  APLCON b("Fit B");
  Vec vec1b = vec1a;
  Vec vec2b = vec2a;
  Vec vec3b = vec3a;

  // for instance b, we link all 4 components at once
  constexpr auto linker4   = [] (Vec& v) -> vector<double*> { return {&v.E, &v.px, &v.py, &v.pz}; };
  b.LinkVariable("Vec1", linker4(vec1b), sigma1);
  b.LinkVariable("Vec2", linker4(vec2b), sigma2);
  b.LinkVariable("Vec3", linker4(vec3b), sigma3);



  // covariances can contain NaN to indicate that they should be kept 0 (so no correlation)
  // can also use APLCON::NaN, which is the identical expression but easier to remember
  constexpr double NaN = numeric_limits<double>::quiet_NaN();

  // we show here how to link covariances.
  // In case you want to set those linked covariances to 0, provide nullptr instead of NaN

  // make some space for a linked covariance between Vec1 components
  // set it to the same values as for instance a
  vector<double> linked_covariance_vec1 = {

  };

  auto equal_vector = [] (const vector<double>& a, const vector<double>& b) -> vector<double> {
    // one should check if sizes of a and b are equal,
    // again that's something APLCON can't do for you
    vector<double> r(a.size());
    for(size_t i=0;i<a.size();i++)
      r[i] = a[i]-b[i];
    return r;
  };

  // finally, do the fit
  // note that many setup exceptions are only thrown here,
  // because only with a fully setup instance it's possible to check many things
  cout.precision(3); // set precision globally, which makes output nicer
  cout << a.DoFit() << endl;

  cout << "Please note that the above fit result might not be meaningful due to totally guessed covariances." << endl;
}
