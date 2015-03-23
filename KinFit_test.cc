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

  k.AddConstraint("BF_e_equal", [] (const map<string, double>& X) {
    return X.at("BF_e_A") - X.at("BF_e_B");
  });
  k.AddConstraint("BF_tau_equal", [] (const map<string, double>& X) {
    return X.at("BF_tau_A") - X.at("BF_tau_B");
  });

}
