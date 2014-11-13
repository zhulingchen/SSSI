// Eikonal solver in C code

#include <math.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"

// function T = eikonal2d(V, dx, sz, sx, iMax)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize nz = mxGetM(prhs[0]);
    mwSize nx = mxGetN(prhs[0]);
    mexPrintf("nz = %d, nx = %d\n", nz, nx);
    
    mxArray *mSlow = mxCreateDoubleMatrix(nz, nx, mxREAL);
    double *pV = mxGetPr(prhs[0]);
    double *pSlow = mxGetPr(mSlow);
    
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            pSlow[i*nx+j] = 1 / pV[i*nx+j];
            mexPrintf("%f\n", pV[i*nx+j]);
        }   
    }
            
    
}

