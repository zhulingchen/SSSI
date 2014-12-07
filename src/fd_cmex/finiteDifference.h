#ifndef _FINITEDIFFERENCE_H
#define _FINITEDIFFERENCE_H

/* ======================================================================
 *
 * dCoef
 * Calculates the coefficients used for finite difference method
 *
 ====================================================================== */
mxArray* dCoef(int order, const char* type);

/* ======================================================================
 *
 * diffOperator
 * Performs higher-order approximation of staggered-grid finite difference
 *
 ====================================================================== */
mxArray* diffOperator(const mxArray *data, const mxArray *coeff, double dist, int dim);




#endif