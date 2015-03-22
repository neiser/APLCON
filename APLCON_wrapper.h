#ifndef __APLCON_wrapper_h__
#define __APLCON_wrapper_h__

// main routines
void c_aplcon_aplcon(const int NVAR, const int MCST);
void c_aplcon_aploop(double X[], double VX[], double F[], int* IRET);
void c_aplcon_chndpv(float* CHI2, int* ND, float* pval);
void c_aplcon_apstat(double* FOPT, int* NFUN, int* NITER);
void c_aplcon_simsel(const double X[], const double VX[], const int NY, const int LIST[], double Y[], double VY[]);
void c_aplcon_simtrn(double X[], double VX[], const int NX);

// printout related routines (useful for debugging)
void c_aplcon_apname(const int I, const char* NAME);
void c_aplcon_ciprv(const int LUP, const double X[], const double VX[], const int N);
void c_aplcon_cfprv(const int LUP, const double X[], const double VX[], const int N);
void c_aplcon_cfcorr(const int LUP, const double V[], const int N);

// setup/config routines
void c_aplcon_aprint(const int LUNP, const int IPR);
void c_aplcon_apdeps(const double EPSF);
void c_aplcon_apderf(const double DERFAC);
void c_aplcon_apderu(const double DERUFC);
void c_aplcon_apdlow(const double DERLOW);
void c_aplcon_apiter(const int ITERMX);
void c_aplcon_aprofl(const int I1, const int I2);
void c_aplcon_apstep(const int I, const double STEP);
void c_aplcon_apfix(const int I);
void c_aplcon_aplimt(const int I, const double XLOW, const double XHIG);
void c_aplcon_aptrin(const int I);
void c_aplcon_apoiss(const int I);
void c_aplcon_abinom(const int I);
void c_aplcon_aplogn(const int I);
void c_aplcon_apsqrt(const int I);
void c_aplcon_apower(const int I);
void c_aplcon_aposit(const int I);

#endif
