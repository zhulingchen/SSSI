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

/* ======================================================================
 *
 * diffOperator
 * DAMPPML Generate the model for damping parameter
 * u = x or z, representing the distance between current position (in PML)
 * and PML inner boundary
 *
 ====================================================================== */
mxArray* dampPml(const mxArray *u, const mxArray *v, double L);


#endif