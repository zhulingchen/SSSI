#ifndef _FINITEDIFFERENCE_H
#define _FINITEDIFFERENCE_H

/* ======================================================================
 *
 * dCoef
 * Calculates the coefficients used for finite difference method
 *
 ====================================================================== */
double* dCoef(int order, const char* type);

/* ======================================================================
 *
 * diffOperator
 * Performs higher-order approximation of staggered-grid finite difference
 *
 ====================================================================== */
/* 2-D case */
double* diffOperator2d(const double *pData, mwSize m, mwSize n, const double *pCoeff, int order, double dist, int dim);
/* 3-D case */
double* diffOperator3d(const double *pData, mwSize n1, mwSize n2, mwSize n3, const double *pCoeff, int order, double dist, int dim);

/* ======================================================================
 *
 * dampPml
 * DAMPPML Generate the model for damping parameter
 * u = x or z, representing the distance between current position (in PML)
 * and PML inner boundary
 *
 ====================================================================== */
double* dampPml(const double *pu, const double *pv, int m, int n, double L);


#endif