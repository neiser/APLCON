#include <APLCON_wrapper.h>
#include <iostream>
#include <vector>

using namespace std;

int main() {
	vector<double> X = {0.1050, 0.1350, 0.0950, 0.1400};
	vector<double> V = {1.0e-4, 0.0e-4, 9.0e-4, 0.0e-4, 0.0e-4, 9.0e-4, 0.0e-4,  0.0e-4,  0.0e-4, 9.0e-4};

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

}
