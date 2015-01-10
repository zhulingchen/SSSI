/* ======================================================================
 *
 * finiteDifference.c
 *
 * Collects all functions used for finite difference method
 *
 * This C source file is free for use in academic research.
 * All rights reserved.
 *
 *
 * Written by Lingchen Zhu (zhulingchen@gmail.com)
 * Center for Signal and Information Processing, Center for Energy & Geo Processing
 * Georgia Institute of Technology
 *
 * ====================================================================== */

#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <math.h>


/* ====================================================================== */
double* dCoef(int order, const char* type)
{
    /* begin of declaration */
    double *pA, *pb, *plhs_mldivide, *pCoeff;
    mxArray *lhs_mldivide[1], *rhs_mldivide[2];
    
    int i, j;
    /* end of declaration */
    
    pA = (double*)mxCalloc(order * order, sizeof(double));
    pb = (double*)mxCalloc(order, sizeof(double));
    
    if (!strcmp(type, "r"))
    {
        pb[0] = 1.0/2.0;
        for (j = 0; j < order; j++)
            for (i = 0; i < order; i++)
                pA[j * order + i] = pow(j+1, 2*(i+1)-1);
    }
    else if (!strcmp(type, "s"))
    {
        pb[0] = 1;
        for (j = 0; j < order; j++)
            for (i = 0; i < order; i++)
                pA[j * order + i] = pow(2*(j+1)-1, 2*(i+1)-1);
    }
    else
    {
        mexErrMsgTxt("Type must be \'s\' or \'r\'!");
    }
    
    /* c = A \ b; */
    /* Point mxArray* to double* */
    rhs_mldivide[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetData(rhs_mldivide[0], pA);
    mxSetM(rhs_mldivide[0], order);
    mxSetN(rhs_mldivide[0], order);
    
    rhs_mldivide[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetData(rhs_mldivide[1], pb);
    mxSetM(rhs_mldivide[1], order);
    mxSetN(rhs_mldivide[1], 1);
    
    /* Call Matlab function to calculate c = A \ b, or equivalently, c = mldivide(A, b) */
    mexCallMATLAB(1, lhs_mldivide, 2, rhs_mldivide, "mldivide");
    
    /* get data for output from mxArray */
    plhs_mldivide = mxGetPr(lhs_mldivide[0]);
    pCoeff = (double*)mxCalloc(order, sizeof(double));
    for (i = 0; i < order; i++)
        pCoeff[i] = plhs_mldivide[i];
    
    /* ATTENTION: Don't forget to free dynamic memory allocated by the mxCreate* function(s) (except for output arrays), otherwise memory leak will occur */
    mxDestroyArray(rhs_mldivide[0]);
    mxDestroyArray(rhs_mldivide[1]);
    mxDestroyArray(lhs_mldivide[0]);
    
    return pCoeff;
}


/* ====================================================================== */
/* 2-D case */
double* diffOperator2d(const double *pData, mwSize m, mwSize n, const double *pCoeff, int order, double dist, int dim)
{
    /* begin of declaration */
    double *pDiffData;
    
    int l, iOrder, i, j, k, idx1, idx2;
    /* end of declaration */
    
    if (dim < 1)
        dim = 1;
    
    l = 2 * order - 1;
    
    if (dim == 1)
    {
        pDiffData = (double*)mxCalloc((m-l) * n, sizeof(double));
        for (iOrder = 0; iOrder < order; iOrder++)
        {
            /*
             * idx1 = ((order+1):(m-l+order)) + (ii-1);
             * idx2 = ((order+1):(m-l+order)) - ii;
             * diffData = diffData + c(ii) * (data(idx1, :) - data(idx2, :)) / d;
             */
            for (j = 0; j < n; j++)
            {
                i = 0;
                idx2 = order - (iOrder + 1);
                for (idx1 = order+iOrder; idx1 <= m-order+iOrder; idx1++)
                {
                    pDiffData[j*(m-l) + i] = pDiffData[j*(m-l) + i] +
                            pCoeff[iOrder] * (pData[j*m + idx1] - pData[j*m + idx2]) / dist;
                    i++;
                    idx2++;
                }
            }
        }
    }
    else    /* dim == 2 */
    {
        pDiffData = (double*)mxCalloc(m * (n-l), sizeof(double));
        for (iOrder = 0; iOrder < order; iOrder++)
        {
            /*
             * idx1 = ((order+1):(n-l+order)) + (ii-1);
             * idx2 = ((order+1):(n-l+order)) - ii;
             * diffData = diffData + c(ii) * (data(:, idx1) - data(:, idx2)) / d;
             */
            for (i = 0; i < m; i++)
            {
                j = 0;
                idx2 = order - (iOrder + 1);
                for (idx1 = order+iOrder; idx1 <= n-order+iOrder; idx1++)
                {
                    pDiffData[j*m + i] = pDiffData[j*m + i] +
                            pCoeff[iOrder] * (pData[idx1*m + i] - pData[idx2*m + i]) / dist;
                    j++;
                    idx2++;
                }
            }
        }
    }
    
    return pDiffData;
}

/* 3-D case */
double* diffOperator3d(const double *pData, mwSize n1, mwSize n2, mwSize n3, const double *pCoeff, int order, double dist, int dim)
{
    /* begin of declaration */
    double *pDiffData;
    
    int l, iOrder, i, j, k, idx1, idx2;
    /* end of declaration */
    
    if (dim < 1)
        dim = 1;
    
    l = 2 * order - 1;
    
    if (dim == 1)
    {
        pDiffData = (double*)mxCalloc((n1-l) * n2 * n3, sizeof(double));
        for (iOrder = 0; iOrder < order; iOrder++)
        {
            /*
             * idx1 = ((order+1):(n1-l+order)) + (ii-1);
             * idx2 = ((order+1):(n1-l+order)) - ii;
             * diffData = diffData + c(ii) * (data(idx1, :, :) - data(idx2, :, :)) / d;
             */
            for (k = 0; k < n3; k++)
            {
                for (j = 0; j < n2; j++)
                {
                    i = 0;
                    idx2 = order - (iOrder + 1);
                    for (idx1 = order+iOrder; idx1 <= n1-order+iOrder; idx1++)
                    {
                        pDiffData[k*((n1-l)*n2) + j*(n1-l) + i] = pDiffData[k*((n1-l)*n2) + j*(n1-l) + i] +
                                pCoeff[iOrder] * (pData[k*(n1*n2) + j*n1 + idx1] - pData[k*(n1*n2) + j*n1 + idx2]) / dist;
                        i++;
                        idx2++;
                    }
                }
            }
        }
    }
    else if (dim == 2)
    {
        pDiffData = (double*)mxCalloc(n1 * (n2-l) * n3, sizeof(double));
        for (iOrder = 0; iOrder < order; iOrder++)
        {
            /*
             * idx1 = ((order+1):(n2-l+order)) + (ii-1);
             * idx2 = ((order+1):(n2-l+order)) - ii;
             * diffData = diffData + c(ii) * (data(:, idx1, :) - data(:, idx2, :)) / d;
             */
            for (k = 0; k < n3; k++)
            {
                for (i = 0; i < n1; i++)
                {
                    j = 0;
                    idx2 = order - (iOrder + 1);
                    for (idx1 = order+iOrder; idx1 <= n2-order+iOrder; idx1++)
                    {
                        pDiffData[k*(n1*(n2-l)) + j*n1 + i] = pDiffData[k*(n1*(n2-l)) + j*n1 + i] +
                                pCoeff[iOrder] * (pData[k*(n1*n2) + idx1*n1 + i] - pData[k*(n1*n2) + idx2*n1 + i]) / dist;
                        j++;
                        idx2++;
                    }
                }
            }
        }
    }
    else    /* dim == 3 */
    {
        pDiffData = (double*)mxCalloc(n1 * n2 * (n3-l), sizeof(double));
        for (iOrder = 0; iOrder < order; iOrder++)
        {
            /*
             * idx1 = ((order+1):(n3-l+order)) + (ii-1);
             * idx2 = ((order+1):(n3-l+order)) - ii;
             * diffData = diffData + c(ii) * (data(:, :, idx1) - data(:, :, idx2)) / d;
             */
            for (j = 0; j < n2; j++)
            {
                for (i = 0; i < n1; i++)
                {
                    k = 0;
                    idx2 = order - (iOrder + 1);
                    for (idx1 = order+iOrder; idx1 <= n3-order+iOrder; idx1++)
                    {
                        pDiffData[k*(n1*n2) + j*n1 + i] = pDiffData[k*(n1*n2) + j*n1 + i] +
                                pCoeff[iOrder] * (pData[idx1*(n1*n2) + j*n1 + i] - pData[idx2*(n1*n2) + j*n1 + i]) / dist;
                        k++;
                        idx2++;
                    }
                }
            }
        }
    }
    
    return pDiffData;
}


/* ====================================================================== */
double* dampPml(const double *pu, const double *pv, mwSize m, mwSize n, double L)
{
    /* begin of declaration */
    double *pd0, *pd;
    
    int i, j;
    /* end of declaration */
    
    const double R = 1e-6;
    const double logR = log(R);
    
    /* d0 = -(3 * v)/(2 * L) * log(R); */
    pd0 = (double*)mxCalloc(m * n, sizeof(double));
    for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
            pd0[j * m + i] = -(3.0 * pv[j * m + i])/(2 * L) * logR;
    
    /* d = d0 .* (u / L).^2; */
    pd = (double*)mxCalloc(m * n, sizeof(double));
    for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
            pd[j * m + i] = pd0[j * m + i] * (pu[j * m + i] / L) * (pu[j * m + i] / L);
    
    /* ATTENTION: Don't forget to free dynamic memory allocated by the mxCalloc function (except for output arrays), otherwise memory leak will occur */
    mxFree(pd0);
    
    return pd;
}