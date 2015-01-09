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
 ====================================================================== */

#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <math.h>


/* ====================================================================== */
double* dCoef(int order, const char* type)
{
    /* begin of declaration */    
    double *pA, *pb, *pCoeff;
    mxArray *lhs_mldivide[1], *rhs_mldivide[2];
    
    int i, j;
    /* end of declaration */
    
    pA = (double*)mxCalloc(order * order, sizeof(double));
    pb = (double*)mxCalloc(order, sizeof(double));
    
    if (!strcmp(type, "r"))
    {
        pb[0] = 1.0/2.0;
        for (j = 0; j < order; j++)
        {
            for (i = 0; i < order; i++)
            {
                pA[j * order + i] = pow(j+1, 2*(i+1)-1);
            }
        }
    }
    else if (!strcmp(type, "s"))
    {
        pb[0] = 1;
        for (j = 0; j < order; j++)
        {
            for (i = 0; i < order; i++)
            {
                pA[j * order + i] = pow(2*(j+1)-1, 2*(i+1)-1);
            }
        }
    }
    else
    {
        mexErrMsgTxt("Type must be \'s\' or \'r\'!");
    }
    
    /* c = A \ b; */
    /* Point mxArray to dynamicData */
    mxSetData(rhs_mldivide[0], pA);
    mxSetM(rhs_mldivide[0], order);
    mxSetN(rhs_mldivide[0], order);
    
    mxSetData(rhs_mldivide[1], pb);
    mxSetM(rhs_mldivide[1], order);
    mxSetN(rhs_mldivide[1], 1);
    
    mexCallMATLAB(1, lhs_mldivide, 2, rhs_mldivide, "mldivide");
    
    /* ATTENTION: Don't forget to free dynamic memory allocated by the mxCalloc function (except for output arrays), otherwise memory leak will occur */
    mxFree(pA);
    mxFree(pb);
    
    /* get data for output from mxArray */
    pCoeff = (double*)mxGetData(lhs_mldivide[0]);
    return pCoeff;
}


/* ====================================================================== */
double* diffOperator(const double *pData, int ndims, const mwSize *pDimsData, const double *pCoeff, int order, double dist, int dim)
{
    /* begin of declaration */
    double *pDiffData;
    
    int l, iOrder, i, j, k, idx1, idx2;

    mwSize n1, n2, n3;
    /* end of declaration */
    
    if (dim < 1)
        dim = 1;
    
    l = 2 * order - 1;
    
    if (ndims <= 2)     /* 2-D case */
    {
        n1 = pDimsData[0];
        n2 = pDimsData[1];
        if (dim == 1)
        {
            pDiffData = (double*)mxCalloc((n1-l) * n2, sizeof(double));
            for (iOrder = 0; iOrder < order; iOrder++)
            {
                /*
                 * idx1 = ((order+1):(n1-l+order)) + (ii-1);
                 * idx2 = ((order+1):(n1-l+order)) - ii;
                 * diffData = diffData + c(ii) * (data(idx1, :) - data(idx2, :)) / d;
                 */
                for (j = 0; j < n2; j++)
                {
                    i = 0;
                    idx2 = order - (iOrder + 1);
                    for (idx1 = order+iOrder; idx1 <= n1-order+iOrder; idx1++)
                    {
                        pDiffData[j*(n1-l) + i] = pDiffData[j*(n1-l) + i] +
                                pCoeff[iOrder] * (pData[j*n1 + idx1] - pData[j*n1 + idx2]) / dist;
                        i++;
                        idx2++;
                    }
                }
            }
        }
        else    /* dim == 2 */
        {
            pDiffData = (double*)mxCalloc(n1 * (n2-l), sizeof(double));
            for (iOrder = 0; iOrder < order; iOrder++)
            {
                /*
                 * idx1 = ((order+1):(n2-l+order)) + (ii-1);
                 * idx2 = ((order+1):(n2-l+order)) - ii;
                 * diffData = diffData + c(ii) * (data(:, idx1) - data(:, idx2)) / d;
                 */
                for (i = 0; i < n1; i++)
                {
                    j = 0;
                    idx2 = order - (iOrder + 1);
                    for (idx1 = order+iOrder; idx1 <= n2-order+iOrder; idx1++)
                    {
                        pDiffData[j*n1 + i] = pDiffData[j*n1 + i] +
                                pCoeff[iOrder] * (pData[idx1*n1 + i] - pData[idx2*n1 + i]) / dist;
                        j++;
                        idx2++;
                    }
                }
            }
        }
    }
    else                /* 3-D case */
    {
        n1 = pDimsData[0];
        n2 = pDimsData[1];
        n3 = pDimsData[2];
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