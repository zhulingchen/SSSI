/* ======================================================================
 *
 * fwdTimeCpmlFor2dAw_mex.c
 *
 * Simulates 2-d acoustic wave forward propagation using finite difference
 * in time domain with partial differential equation (PDE)
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
#include "finiteDifference.h"

// function [data, snapshot] = fwdTimeCpmlFor2dAw(v, source, nDiffOrder, nBoundary, dz, dx, dt)

// input arguments
#define VM_IN           prhs[0]
#define SOURCE_IN       prhs[1]
#define DIFFORDER_IN	prhs[2]
#define BOUNDARY_IN     prhs[3]
#define DZ_IN           prhs[4]
#define DX_IN           prhs[5]
#define DT_IN           prhs[6]

// output arguments
#define DATA_OUT        plhs[0]
#define SNAPSHOT_OUT    plhs[1]

// the gateway routine
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *pVelocityModel, *pSource, *pData, *pSnapshot, *pCoeff, *pOldFdm, *pCurFdm, *pNewFdm;
    int diffOrder, boundary, dz, dx, dt, nz, nx, nt, l;
    const mwSize *pSourceDims;
    
    mxArray *coeff, *oldFdm, *curFdm, *newFdm;
    
    if (nrhs < 7)
        mexErrMsgTxt("All 7 input arguments shall be provided!");
    
    // ATTENTION: mxGetPr might just produce a 1D array that is linearized according to Matlab convention (column order)
    pVelocityModel = mxGetPr(VM_IN);
    pSource = mxGetPr(SOURCE_IN);
    diffOrder = *mxGetPr(DIFFORDER_IN);
    boundary = *mxGetPr(BOUNDARY_IN);
    dz = *mxGetPr(DZ_IN);
    dx = *mxGetPr(DX_IN);
    dt = *mxGetPr(DT_IN);
    
    pSourceDims = mxGetDimensions(SOURCE_IN);
    nz = pSourceDims[0];
    nx = pSourceDims[1];
    nt = pSourceDims[2];
    
    //mexPrintf("nz = %d, nx = %d, nt = %d\n", nz, nx, nt);
    DATA_OUT = mxCreateDoubleMatrix(nx, nt, mxREAL);
    pData = mxGetPr(DATA_OUT);
    
    mwSize pDimsSnapshot[3] = {nz, nx, nt};
    SNAPSHOT_OUT = mxCreateNumericArray(3, pDimsSnapshot, mxDOUBLE_CLASS, mxREAL);
    pSnapshot = mxGetPr(SNAPSHOT_OUT);
    
    coeff = dCoef(diffOrder, "s");
    pCoeff = mxGetPr(coeff);
    l = 2 * diffOrder - 1;
    oldFdm = mxCreateDoubleMatrix(nz+2*l, nx+2*l, mxREAL);
    curFdm = mxCreateDoubleMatrix(nz+2*l, nx+2*l, mxREAL);
    newFdm = mxCreateDoubleMatrix(nz+2*l, nx+2*l, mxREAL);
    
    
    // ATTENTION: Don't forget to free dynamic memory allocated by MXCREATE* functions (except for output arrays), otherwise memory leak will occur
    mxDestroyArray(coeff);
    mxDestroyArray(oldFdm);
    mxDestroyArray(curFdm);
    mxDestroyArray(newFdm);
}