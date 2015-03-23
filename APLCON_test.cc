#include <iostream>
#include <vector>

extern "C" {
#include <APLCON_wrapper.h>
}

using namespace std;

int main() {
	vector<double> X = {0.1050, 0.1350, 0.0950, 0.1400};
	vector<double> V = {1.0e-4, 0.0e-4, 9.0e-4, 0.0e-4, 0.0e-4, 9.0e-4, 0.0e-4,  0.0e-4,  0.0e-4, 9.0e-4};
	vector<double> F(2);

	int NVAR = 4;
	int NEQS = 2;
	c_aplcon_aplcon(NVAR,NEQS);
	int LUP = 6;
	int IPR = 6;
	c_aplcon_aprint(LUP,IPR);

	c_aplcon_ciprv(LUP,X.data(),V.data(),NVAR); //print initial status
	c_aplcon_cfcorr(LUP,V.data(),NVAR);
	c_aplcon_aprofl(4,0);

	int IRET = -1;
	do {
		F[0] = X[0] - X[1];
		F[1] = X[2] - X[3];
		c_aplcon_aploop(X.data(), V.data(), F.data(), &IRET);
	}
	while(IRET<0);

	float chi2, pval;
	double fopt;
	int ndof, nfun, niter;
  vector<double> pulls(X.size());
	c_aplcon_chndpv(&chi2,&ndof,&pval);
	c_aplcon_apstat(&fopt,&nfun,&niter);
  c_aplcon_appull(pulls.data());
  
	cout <<  "Chi2 = " << chi2
	     << " NDF = " << ndof
	     << " pval = " << pval
	     << " Fopt = " << fopt
	     << " nFun = " << nfun
	     << " nIter = " << niter
	     << endl;
  
  for(auto p : pulls) {
    cout << p << endl;
  }

}
