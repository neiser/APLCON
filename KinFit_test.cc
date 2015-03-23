#include <iostream>
#include <vector>
#include <KinFit.h>

using namespace std;

int main() {
  KinFit k;
  k.AddMeasuredVariable("BF_e_A", 0.1050, 0.01);
  k.AddMeasuredVariable("BF_e_B", 0.135, 0.03);
  k.AddMeasuredVariable("BF_tau_A", 0.095, 0.03);
  k.AddMeasuredVariable("BF_tau_B", 0.14, 0.03);

  /*k.AddConstraint("BF_e_equal", [] (const map<string, double>& X) {
    return X.at("BF_e_A") - X.at("BF_e_B");
  });
  k.AddConstraint("BF_tau_equal", [] (const map<string, double>& X) {
    return X.at("BF_tau_A") - X.at("BF_tau_B");
  });

  k.AddConstraint("BF_e_A", "BF_e_B", [] (double a, double b) {
    return a - b;
  });

  k.AddConstraint("BF_e_A", "BF_e_B", "BF_e_C", [] (double a, double b, double c) {
    return a - b + c;
  });*/

  //const vector<const string> t = ;

  k.AddConstraint("test", {"bla", "blu"},
                  [] (double a, double b) { return a - b;});

  //k.DoFit();
}
