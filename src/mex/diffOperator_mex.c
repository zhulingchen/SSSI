/* ======================================================================
 *
 * diffOperator_mex.c
 *
 * Performs higher-order approximation of staggered-grid finite difference
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

// input arguments
#define DATA_IN     prhs[0]
#define COEFF_IN    prhs[1]
#define DIST_IN     prhs[2]
#define DIM_IN      prhs[3]

// output arguments
#define DIFF_OUT    plhs[0]

// the gateway routine
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *pData, *pCoeff, dist, *pDiffData;
    int dim, l, n1, n2, n3;
    mwSize ndims, order;
    const mwSize *pDims;
    
    int iOrder, i, j, k, idx1, idx2;
    
    if (nrhs < 3)
    {
        mexErrMsgTxt("Input data, finite difference coefficients and distance per sample should be all provided!");
    }
    if (nrhs < 4)
    {
        dim = 1;
    }
    
    // ATTENTION: mxGetPr might just produce a 1D array that is linearized according to Matlab convention (column order)
    pData = mxGetPr(DATA_IN);
    pCoeff = mxGetPr(COEFF_IN);
    dist = *mxGetPr(DIST_IN);
    dim = *mxGetPr(DIM_IN);
    
    ndims = mxGetNumberOfDimensions(DATA_IN);
    order = mxGetNumberOfElements(COEFF_IN);
    l = 2 * order - 1;
    
    pDims = mxGetDimensions(DATA_IN);
    
    if (ndims <= 2)     // 2-D case
    {
        n1 = pDims[0];
        n2 = pDims[1];
        if (dim == 1)
        {
            DIFF_OUT = mxCreateDoubleMatrix(n1-l, n2, mxREAL);
            pDiffData = mxGetPr(DIFF_OUT);
            for (iOrder = 0; iOrder < order; iOrder++)
            {
                //idx1 = ((order+1):(n1-l+order)) + (ii-1);
                //idx2 = ((order+1):(n1-l+order)) - ii;
                //diffData = diffData + c(ii) * (data(idx1, :) - data(idx2, :)) / d;
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
        else    // dim == 2
        {
            DIFF_OUT = mxCreateDoubleMatrix(n1, n2-l, mxREAL);
            pDiffData = mxGetPr(DIFF_OUT);
            for (iOrder = 0; iOrder < order; iOrder++)
            {
                //idx1 = ((order+1):(n2-l+order)) + (ii-1);
                //idx2 = ((order+1):(n2-l+order)) - ii;
                //diffData = diffData + c(ii) * (data(:, idx1) - data(:, idx2)) / d;
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
    else                // 3-D case
    {
        n1 = pDims[0];
        n2 = pDims[1];
        n3 = pDims[2];
        if (dim == 1)
        {
            mwSize pDimsNew[3] = {n1-l, n2, n3};
            DIFF_OUT = mxCreateNumericArray(ndims, pDimsNew, mxDOUBLE_CLASS, mxREAL);
            pDiffData = mxGetPr(DIFF_OUT);
            for (iOrder = 0; iOrder < order; iOrder++)
            {
                //idx1 = ((order+1):(n1-l+order)) + (ii-1);
                //idx2 = ((order+1):(n1-l+order)) - ii;
                //diffData = diffData + c(ii) * (data(idx1, :, :) - data(idx2, :, :)) / d;
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
            mwSize pDimsNew[3] = {n1, n2-l, n3};
            DIFF_OUT = mxCreateNumericArray(ndims, pDimsNew, mxDOUBLE_CLASS, mxREAL);
            pDiffData = mxGetPr(DIFF_OUT);
            for (iOrder = 0; iOrder < order; iOrder++)
            {
                //idx1 = ((order+1):(n2-l+order)) + (ii-1);
                //idx2 = ((order+1):(n2-l+order)) - ii;
                //diffData = diffData + c(ii) * (data(:, idx1, :) - data(:, idx2, :)) / d;
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
        else    // dim == 3
        {
            mwSize pDimsNew[3] = {n1, n2, n3-l};
            DIFF_OUT = mxCreateNumericArray(ndims, pDimsNew, mxDOUBLE_CLASS, mxREAL);
            pDiffData = mxGetPr(DIFF_OUT);
            for (iOrder = 0; iOrder < order; iOrder++)
            {
                //idx1 = ((order+1):(n3-l+order)) + (ii-1);
                //idx2 = ((order+1):(n3-l+order)) - ii;
                //diffData = diffData + c(ii) * (data(:, :, idx1) - data(:, :, idx2)) / d;
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
}