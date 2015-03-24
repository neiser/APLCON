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

  auto make_equal = [] (double a, double b) { return a - b;};
  k.AddConstraint("BF_e_equal", {"BF_e_A", "BF_e_B"}, make_equal);
  k.AddConstraint("BF_tau_equal", {"BF_tau_A", "BF_tau_B"}, make_equal);

  k.DoFit();
}
