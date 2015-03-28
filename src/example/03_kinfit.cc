#include <iostream>
#include <APLCON.hpp>

using namespace std;

int main() {


  APLCON a("KinFit");

  struct Vec {
    double E;
    double px;
    double py;
    double pz;
  };

  Vec vec1 = { 1, 2, 3, 4};
  Vec vec2 = { 5, 6, 7, 8};
  Vec vec3 = { 9,10,11,12};
  Vec vec4 = {13,14,15,16};

  auto linker   = [] (Vec& v) -> vector<double*> { return {&v.E, &v.px, &v.py, &v.pz}; };
  auto linker_E = [] (Vec& v) -> vector<double*> { return {&v.E}; };
  auto linker_p = [] (Vec& v) -> vector<double*> { return {&v.px, &v.py, &v.pz}; };

  a.LinkVariable("Vec1_E", linker_E(vec1), vector<double>{1});
  a.LinkVariable("Vec1_p", linker_p(vec1), vector<double>{2});
  a.LinkVariable("Vec2_E", linker_E(vec2), vector<double>{3});
  a.LinkVariable("Vec2_p", linker_p(vec2), vector<double>{4});
  a.LinkVariable("Vec3",   linker(vec3),   vector<double>{5});
  a.LinkVariable("Vec4",   linker(vec4),   vector<double>{6});


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


//  auto scalar2scalar = [] (const double& a) -> double {
//    return  a;
//  };

//  auto scalar2vector = [] (const double& a) -> vector<double> {
//    return  {a};
//  };

//  auto vector2scalar = [] (const vector<double>& a) -> double {
//    return  a[0] - a[1];
//  };

//  auto vector2vector = [] (const vector<double>& a) -> vector<double> {
//    return  {a[0] - a[1]};
//  };

//  a.AddConstraint("1", {"BF_e"}, scalar2scalar);
//  a.AddConstraint("2", {"BF_e"}, scalar2vector);
//  a.AddConstraint("3", {"BF_e"}, vector2scalar);
//  a.AddConstraint("4", {"BF_e"}, vector2vector);

  //  auto wrong = [] (const vector<double>& a, double b) -> vector<double> {
  //    return  {a[0] - a[1] + b};
  //  };

    //a.AddConstraint("4", {"BF_e", "BF_tau"}, wrong);
}
