#ifndef __APLCON_wrapper_h__
#define __APLCON_wrapper_h__

extern "C" {
	void c_aplcon_aplcon(int NVAR, int MCST);
	void c_aplcon_aprint(int LUNP, int IPR);
	void c_aplcon_apname(int I, const char* NAME);
	void c_aplcon_ciprv(int LUP, const double X[], const double VX[], int N);  
}

#endif
