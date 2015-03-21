#ifndef __APLCON_wrapper_h__
#define __APLCON_wrapper_h__

extern "C" {
	void c_aplcon_aplcon(int NVAR, int MCST);
	void c_aplcon_apname(int I, const char* NAME);
}

#endif
