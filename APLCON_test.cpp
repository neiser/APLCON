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

	c_aplcon_apname(1,"BR to e (A)");
	c_aplcon_apname(2,"BR to e (B)");
	c_aplcon_apname(3,"BR to tau (A)");
	c_aplcon_apname(4,"BR to tau (B)");
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
}
