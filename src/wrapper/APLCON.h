#ifndef APLCON_wrapper_h
#define APLCON_wrapper_h

// main routines
/**
 * @brief Initialize APLCON
 * @param NVAR number of variables
 * @param MCST number of constraints
 */
void c_aplcon_aplcon(const int NVAR, const int MCST);
/**
 * @brief Loop function for fit
 * @param X current variable values
 * @param VX current covariance matrix (symmetrized)
 * @param F current constraint values
 * @param IRET status of fit iteration
 */
void c_aplcon_aploop(double X[], double VX[], double F[], int* IRET);

// routines to obtain results
/**
 * @brief Obtain ChiSquare, NDoF, Probability
 * @see c_aplcon_apstat
 * @param CHI2 ChiSquare as float
 * @param ND Number of degrees of freedom
 * @param PVAL Probability
 */
void c_aplcon_chndpv(float* CHI2, int* ND, float* PVAL);
/**
 * @brief Obtain information after fit
 * @param FOPT ChiSqare as double
 * @param NFUN number of function calls
 * @param NITER number of iterations
 */
void c_aplcon_apstat(double* FOPT, int* NFUN, int* NITER);
/**
 * @brief Obtain pulls
 * @param PULLS Array of pulls for each variable in X
 */
void c_aplcon_appull(double* PULLS);

// variable reduction (currently unused)
//void c_aplcon_simsel(const double X[], const double VX[], const int NY, const int LIST[], double Y[], double VY[]);
//void c_aplcon_simtrn(double X[], double VX[], const int NX);

// printout related routines (useful for debugging, but unused now)
//void c_aplcon_ciprv(const int LUP, const double X[], const double VX[], const int N);
//void c_aplcon_cfprv(const int LUP, const double X[], const double VX[], const int N);
//void c_aplcon_cfcorr(const int LUP, const double V[], const int N);

// setup/config routines
/**
 * @brief Setup verbosity
 * @param LUNP device to printout messages
 * @param IPR verbosity level
 */
void c_aplcon_aprint(const int LUNP, const int IPR);
/**
 * @brief Setup constraint accuracy
 * @param EPSF Constraint accuracy
 */
void c_aplcon_apdeps(const double EPSF);
/**
 * @brief Setup chi2 accuracy
 * @param EPSCHI chi2 accuracy
 */
void c_aplcon_apepschi(const double EPSCHI);
/**
 * @brief Setup measured stepsize factor
 * @param DERFAC measured stepsize factor
 */
void c_aplcon_apderf(const double DERFAC);
/**
 * @brief Setup unmeasured stepsize factor
 * @param DERUFC unmeasured stepsize factor
 */
void c_aplcon_apderu(const double DERUFC);
/**
 * @brief Setup minimal stepsize factor
 * @param DERLOW minimal stepsize factor
 */
void c_aplcon_apdlow(const double DERLOW);
/**
 * @brief Setup maximum number of iterations
 * @param ITERMX number of maximum iterations
 */
void c_aplcon_apiter(const int ITERMX);
/**
 * @brief Setup stepsize for variable with index I
 * @param I index of variable (starting from 1)
 * @param STEP stepsize of variable, 0 for fixed
 */
void c_aplcon_apstep(const int I, const double STEP);
/**
 * @brief Set variable I to be fixed
 * @param I index of variable
 */
void c_aplcon_apfix(const int I);
/**
 * @brief Set high and low limit of variable I
 * @param I index of variable
 * @param XLOW lower limit
 * @param XHIG upper limit
 */
void c_aplcon_aplimt(const int I, const double XLOW, const double XHIG);
/**
 * @brief Setup variable I to be poissionian distributed
 * @param I index of variable
 */
void c_aplcon_apoiss(const int I);
/**
 * @brief Setup variable to be sqrt-transformed
 * @param I index of variable
 */
void c_aplcon_apsqrt(const int I);
/**
 * @brief Setup variable to be log-normal distributed
 * @param I index of variable
 */
void c_aplcon_aplogn(const int I);

// rather undocumented additional APLCON routines
//void c_aplcon_abinom(const int I);
//void c_aplcon_apower(const int I);
//void c_aplcon_aposit(const int I);
//void c_aplcon_aprofl(const int I1, const int I2);
//void c_aplcon_aptrin(const int I);

#endif
