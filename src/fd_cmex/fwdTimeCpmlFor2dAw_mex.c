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
#include <math.h>
#include <string.h>

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
#define TEST_OUT        plhs[2]

// the gateway routine
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // begin of declaration
    double *pVelocityModel, *pSource, *pData, *pSnapshot;
    double dz, dx, dt;
    int diffOrder, boundary;
    
    int l, i, j;
    mwSize nz, nx, nt;
    const mwSize *pDimsSource;
    mwSize pDimsSnapshot[3] = {0};
    
    mxArray *coeff, *oldFdm, *curFdm, *newFdm;
    double *pCoeff, *pOldFdm, *pCurFdm, *pNewFdm;
    
    mxArray *uDampLeft, *vDampLeft, *uDampRight, *vDampRight, *uDampDown, *vDampDown;
    double *puDampLeft, *pvDampLeft, *puDampRight, *pvDampRight, *puDampDown, *pvDampDown;
    
    mxArray *xDampLeft, *xDampRight, *xDamp, *xb, *zDampDown, *zDamp, *zb;
    double *pxDampLeft, *pxDampRight, *pxDamp, *pxb, *pzDampDown, *pzDamp, *pzb;
    
    mxArray *vdtSq;
    double *pVdtSq;

    mxArray *zPhi, *xPhi, *zA, *xA, *zPsi, *xPsi, *zP, *xP;
    double *pzPhi, *pxPhi, *pzA, *pxA, *pzPsi, *pxPsi, *pzP, *pxP;
    // end of declaration
    
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
    
    pDimsSource = mxGetDimensions(SOURCE_IN);
    nz = pDimsSource[0];
    nx = pDimsSource[1];
    nt = pDimsSource[2];
    mxAssert(nz == mxGetM(VM_IN), "Velocity model and source grids should have the same z-axis grids!");
    mxAssert(nx == mxGetN(VM_IN), "Velocity model and source grids should have the same x-axis grids!");
    
    // initialize storage
    DATA_OUT = mxCreateDoubleMatrix(nx, nt, mxREAL);
    pData = mxGetPr(DATA_OUT);
    
    pDimsSnapshot[0] = nz;
    pDimsSnapshot[1] = nx;
    pDimsSnapshot[2] = nt;
    SNAPSHOT_OUT = mxCreateNumericArray(3, pDimsSnapshot, mxDOUBLE_CLASS, mxREAL);
    pSnapshot = mxGetPr(SNAPSHOT_OUT);
    
    coeff = dCoef(diffOrder, "s");
    pCoeff = mxGetPr(coeff);
    l = 2 * diffOrder - 1;
    oldFdm = mxCreateDoubleMatrix(nz+2*l, nx+2*l, mxREAL);
    curFdm = mxCreateDoubleMatrix(nz+2*l, nx+2*l, mxREAL);
    newFdm = mxCreateDoubleMatrix(nz+2*l, nx+2*l, mxREAL);
    
    // damp profile of x-axis
    uDampLeft = mxCreateDoubleMatrix(nz, boundary, mxREAL);
    puDampLeft = mxGetPr(uDampLeft);
    for (j = 0; j < boundary; j++)
        for (i = 0; i < nz; i++)
            puDampLeft[j * nz + i] = (boundary - j) * dx;
    vDampLeft = mxCreateDoubleMatrix(nz, boundary, mxREAL);
    pvDampLeft = mxGetPr(vDampLeft);
    memcpy(pvDampLeft, pVelocityModel, sizeof(double) * nz * boundary);
    xDampLeft = dampPml(uDampLeft, vDampLeft, boundary * dx);
    pxDampLeft = mxGetPr(xDampLeft);
    
    uDampRight = mxCreateDoubleMatrix(nz, boundary, mxREAL);
    puDampRight = mxGetPr(uDampRight);
    for (j = 0; j < boundary; j++)
        for (i = 0; i < nz; i++)
            puDampRight[j * nz + i] = (j + 1) * dx;
    vDampRight = mxCreateDoubleMatrix(nz, boundary, mxREAL);
    pvDampRight = mxGetPr(vDampRight);
    memcpy(pvDampRight, pVelocityModel + nz * (nx-boundary), sizeof(double) * nz * boundary);
    xDampRight = dampPml(uDampRight, vDampRight, boundary * dx);
    pxDampRight = mxGetPr(xDampRight);
    
    xDamp = mxCreateDoubleMatrix(nz, nx, mxREAL);
    pxDamp = mxGetPr(xDamp);
    memcpy(pxDamp, pxDampLeft, sizeof(double) * nz * boundary);
    memcpy(pxDamp + nz * (nx-boundary), pxDampRight, sizeof(double) * nz * boundary);
    
    xb = mxCreateDoubleMatrix(nz, nx, mxREAL);
    pxb = mxGetPr(xb);
    for (j = 0; j < nx; j++)
        for (i = 0; i < nz; i++)
            pxb[j * nz + i] = exp(-pxDamp[j * nz + i] * dt);
    
    // damp profile of z-axis
    uDampDown = mxCreateDoubleMatrix(boundary, nx, mxREAL);
    puDampDown = mxGetPr(uDampDown);
    for (j = 0; j < nx; j++)
        for(i = 0; i < boundary; i++)
            puDampDown[j * boundary + i] = (i + 1) * dz;
    vDampDown = mxCreateDoubleMatrix(boundary, nx, mxREAL);
    pvDampDown = mxGetPr(vDampDown);
    for (j = 0; j < nx; j++)
        for(i = 0; i < boundary; i++)
            pvDampDown[j * boundary + i] = pVelocityModel[j * nz + (nz - boundary + i)];
    zDampDown = dampPml(uDampDown, vDampDown, boundary * dz);
    pzDampDown = mxGetPr(zDampDown);
    
    zDamp = mxCreateDoubleMatrix(nz, nx, mxREAL);
    pzDamp = mxGetPr(zDamp);
    for (j = 0; j < nx; j++)
        for (i = nz-boundary; i < nz; i++)
            pzDamp[j * nz + i] = pzDampDown[j * boundary + i-(nz-boundary)];
    
    zb = mxCreateDoubleMatrix(nz, nx, mxREAL);
    pzb = mxGetPr(zb);
    for (j = 0; j < nx; j++)
        for (i = 0; i < nz; i++)
            pzb[j * nz + i] = exp(-pzDamp[j * nz + i] * dt);
    
    /* ======================================================================
     * 2-D Acoustic Wave Forward-Time Modeling
     * ====================================================================== */
    // additional arrays for storage intermediate results
    zPhi = mxCreateDoubleMatrix(nz+2*l, nx, mxREAL);
    xPhi = mxCreateDoubleMatrix(nz, nx+2*l, mxREAL);
    zA = mxCreateDoubleMatrix(nz+2*l, nx, mxREAL);
    xA = mxCreateDoubleMatrix(nz, nx+2*l, mxREAL);
    zPsi = mxCreateDoubleMatrix(nz+l, nx, mxREAL);
    xPsi = mxCreateDoubleMatrix(nz, nx+l, mxREAL);
    zP = mxCreateDoubleMatrix(nz+l, nx, mxREAL);
    xP = mxCreateDoubleMatrix(nz, nx+l, mxREAL);
    
    vdtSq = mxCreateDoubleMatrix(nz, nx, mxREAL);
    pVdtSq = mxGetPr(vdtSq);
    for (j = 0; j < nx; j++)
        for (i = 0; i < nz; i++)
            pVdtSq[j * nz + i] = (pVelocityModel[j * nz + i] * dt) * (pVelocityModel[j * nz + i] * dt);
    
    
    // test begin
    TEST_OUT = zb;
    // test end
    
    
    // ATTENTION: Don't forget to free dynamic memory allocated by MXCREATE* functions (except for output arrays), otherwise memory leak will occur
    mxDestroyArray(coeff);
    mxDestroyArray(oldFdm);
    mxDestroyArray(curFdm);
    mxDestroyArray(newFdm);
    mxDestroyArray(uDampLeft);
    mxDestroyArray(vDampLeft);
    mxDestroyArray(xDampLeft);
    mxDestroyArray(uDampRight);
    mxDestroyArray(vDampRight);
    mxDestroyArray(xDampRight);
    mxDestroyArray(xDamp);
    mxDestroyArray(uDampDown);
    mxDestroyArray(vDampDown);
    mxDestroyArray(zDampDown);
    mxDestroyArray(zDamp);
}