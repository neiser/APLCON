#include <iostream>
#include <vector>
#include <KinFit.h>

using namespace std;

int main() {
  KinFit k;
  k.AddMeasuredVariable("bla", 0.4, 0.03);
  k.AddMeasuredVariable("blu", 0.5, 0.04, KinFit::Distribution_t::Poissonian);
  k.AddMeasuredVariable("bli", 0.6, 0.05, KinFit::Distribution_t::LogNormal, 0.1, 0.5);
  k.AddMeasuredVariable("dummy", 0.6, 0);
}
