#include <iostream>
#include <APLCON.hpp>

using namespace std;

int main() {

  // this example illustrates some more advanced
  // usage of the interface,
  // like setting up (possibly linked) covariances
  // and more complex constraint functions

  // please feel free to modify the code here and see
  // if APLCON will throw you a meaningful exception in
  // case you set up something stupid :)

  APLCON a("Covariances");
  APLCON b("Linked Covariances");

  // as an example structure, we use something
  // which looks like a Lorentz vector
  struct Vec {
    double E;
    double px;
    double py;
    double pz;
  };

  // just some particles..
  Vec vec1a = { 1, 2, 3, 4};
  Vec vec2a = { 5, 6, 7, 8};
  Vec vec1b = vec1a;
  Vec vec2b = vec2a;


  // you're totally free in
  // how the fields of your data structure are linked

  // for instance a, we separate E and p
  auto linker_E = [] (Vec& v) -> vector<double*> { return {&v.E}; };
  auto linker_p = [] (Vec& v) -> vector<double*> { return {&v.px, &v.py, &v.pz}; };
  a.LinkVariable("Vec1_E", linker_E(vec1a), vector<double>{1});
  a.LinkVariable("Vec1_p", linker_p(vec1a), vector<double>{1});
  a.LinkVariable("Vec2_E", linker_E(vec2a), vector<double>{2});
  a.LinkVariable("Vec2_p", linker_p(vec2a), {2});

  // for instance b, we link all 4 components at once
  auto linker4   = [] (Vec& v) -> vector<double*> { return {&v.E, &v.px, &v.py, &v.pz}; };
  b.LinkVariable("Vec1", linker4(vec1b), vector<double>{1});
  b.LinkVariable("Vec2", linker4(vec2b), vector<double>{2});


  // SetCovariance was already used for scalar variables in the first example
  // here we show how to setup covariances for vector-valued variables
  // and how to setup

  a.SetCovariance("Vec1_p", "Vec1_p",
                  vector<double>{
                    // the double slashes mark the diagonal elements
                    //
                    1, //
                    2, 3 //
                  });

  a.SetCovariance("Vec1_E", "Vec1_p", // variables give 1 row and 3 columns
                  vector<double>{
                    4, 5, 6
                  });

  a.SetCovariance("Vec2_p", "Vec2_E", // variables give 3 rows and 1 column
                  vector<double>{
                    4,
                    5,
                    6
                  });

  a.SetCovariance("Vec3", "Vec2_p", // variables give 4 rows and 3 columns
                  vector<double>{
                    7,  8,  9,
                   10, 11, 12,
                   13, 14, 15,
                   16, 17, 18
                  });


  auto equal_energy = [] (const vector<double>& a, const vector<double>& b) -> vector<double> {
    return {a[0] - b[0]}; // might return more than one constraint
  };

  auto equal_momentum_3 = [] (const vector<double>& a, const vector<double>& b) -> vector<double> {
    // one should check that the vectors a, b have the appropiate lengths...
    return {
      a[0] - b[0],
      a[1] - b[1],
      a[2] - b[2]
    }; // might return more than one constraint
  };

  auto equal_vector = [] (const vector<double>& a, const vector<double>& b) -> vector<double> {
    // TODO: check check if sizes of a and b are equal
    vector<double> r(a.size());
    for(size_t i=0;i<a.size();i++)
      r[i] = a[i]-b[i];
    return r;
  };

  a.AddConstraint("equal_energy",    {"Vec1_E", "Vec2_E"}, equal_energy);
  a.AddConstraint("equal_momentum",  {"Vec1_p", "Vec2_p"}, equal_momentum_3);
  a.AddConstraint("equal_fourvector",  {"Vec3", "Vec4"}, equal_vector);

  cout << a.DoFit() << endl;

}
