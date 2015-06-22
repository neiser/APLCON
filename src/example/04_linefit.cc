#include <iostream>
#include <APLCON.hpp>
#include <limits>
#include <iomanip>
#include <functional>

using namespace std;

int main() {

  // this example fits a straight line given by
  // f(x) = a + b*x
  // to datapoints (x, y)
  // with errors in y (and later in x)
  
  struct data_t {
    vector<double> x;
    vector<double> sx; // errors in x
    vector<double> y;
    vector<double> sy; // errors in y    
  };
  
  const data_t data{
    {   1,    2,    3,    4},
    { 0.2, 0.23, 0.16, 0.21},  
    { 1.1, 1.95, 2.02, 3.98},
    {0.08, 0.04, 0.11, 0.07} 
  };
    
  // setup first APLCON
  
  APLCON::Fit_Settings_t settings = APLCON::Fit_Settings_t::Default;
  settings.MaxIterations = 500;
  APLCON f1("StraightLineFit", settings);
  
  // since we deal with many variables, we use the linker interface
  // you need to pay attention to references and so on 
  // to get the proper vector of pointers
  auto linker = [] (vector<double>& v) {
    vector<double*> v_p(v.size());
    transform(v.begin(), v.end(), v_p.begin(), addressof<double>);
    return v_p;
  };
  
  // also save some original values for second fit,
  // since the data is modified due to Linker interface to APLCON
  data_t data1 = data;
  
  APLCON::Variable_Settings_t var_settings = APLCON::Variable_Settings_t::Default;
  var_settings.StepSize = 0; // stepsize=0 means fixed variable
  f1.LinkVariable("x", linker(data1.x), vector<double>{0}, {var_settings}); // use x as fixed variables
  f1.LinkVariable("y", linker(data1.y), linker(data1.sy));
  
  // the two parameters are simply unmeasured
  f1.AddUnmeasuredVariable("a");
  f1.AddUnmeasuredVariable("b");
  
  // now the constraint, which requires
  // a + b*x_i - y_i to be zero for each point (x_i, y_i)
  // this is then a contraint for the parameters a and b

  auto residuals = [] (const vector< vector<double> >& arg) {
    // arg is a bit difficult to access, so define aliases
    // see AddConstraint call why this is the correct order...
    const double& a = arg[0][0];
    const double& b = arg[1][0];
    const vector<double>& x = arg[2];    
    const vector<double>& y = arg[3];
    
    // we use again std::transform to calculate the residuals
    vector<double> residuals(y.size());
    transform(x.begin(), x.end(), y.begin(), residuals.begin(), 
              [&a, &b] (const double& x_i, const double& y_i) {
      return a + b*x_i - y_i;
    });
    return residuals;
  };
  f1.AddConstraint("residuals", vector<string>{"a", "b", "x", "y"}, residuals);
  
  // just output everything after doing the fit
  const APLCON::Result_t& r1 = f1.DoFit();
  cout << r1 << endl;
  
  // again, we can cross-check with 
  // the well-known formulas for straight-line fitting
  double S   = 0;
  double Sx  = 0;
  double Sy  = 0;
  double Sxx = 0;
  double Sxy = 0;
  for(size_t i=0;i<data.x.size();i++) {
    const double s = data.sy[i];
    const double s2 = s*s;
    const double x_  = data.x[i];
    const double y_  = data.y[i];
    S   += 1/s2;
    Sx  += x_/s2;
    Sy  += y_/s2;
    Sxx += x_*x_/s2;
    Sxy += x_*y_/s2;
  }
  const double D = S*Sxx-Sx*Sx;
  double a_direct = (Sxx*Sy - Sx*Sxy)/D;
  double b_direct = (S*Sxy-Sx*Sy)/D;
  
  // print comparison
  cout << "+++++++++++++++++" << endl;
  cout << "Direct calculation: a=" << a_direct 
       << " b=" << b_direct 
       << endl;
  cout << "APLCON says:        a=" << r1.Variables.at("a").Value.After 
       << " b=" << r1.Variables.at("b").Value.After 
       << endl;
  cout << "+++++++++++++++++" << endl << endl;
  
  // Now, APLCON becomes really cool when we start introducing x AND y errors
  // which is very complicated when done explicitly, but easy with APLCON
  
  APLCON f2("StraightLineFitWithXYErrors", settings);
  
  // now link x as well with uncertainties (the only difference compared to above)
  data_t data2 = data;  
  f2.LinkVariable("x", linker(data2.x), linker(data2.sx));  
  f2.LinkVariable("y", linker(data2.y), linker(data2.sy));
  
  // the two parameters are simply unmeasured
  f2.AddUnmeasuredVariable("a");
  f2.AddUnmeasuredVariable("b");
  
  // oh, we can use the same constraint, 
  // because APLCON does the nasty error business
  f2.AddConstraint("residuals", vector<string>{"a", "b", "x", "y"}, residuals);
  
  // output everything
  const APLCON::Result_t& r2 = f2.DoFit();
  cout << r2 << endl;
  
  // compare to normal fit
  cout << "+++++++++++++++++" << endl;
  cout << "Direct calculation:        a=" << a_direct 
       << " b=" << b_direct 
       << endl;
  cout << "APLCON says with x-errors: a=" << r2.Variables.at("a").Value.After 
       << " b=" << r2.Variables.at("b").Value.After 
       << endl;
  cout << "+++++++++++++++++" << endl << endl;
}
