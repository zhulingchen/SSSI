/* ======================================================================
 *
 * rvsTimeCpmlFor2dAw_mex.c
 *
 * Simulates 2-d acoustic wave reverse propagation using finite difference
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

// input arguments
#define VM_IN           prhs[0]
#define DATA_IN         prhs[1]
#define DIFFORDER_IN	prhs[2]
#define BOUNDARY_IN     prhs[3]
#define DZ_IN           prhs[4]
#define DX_IN           prhs[5]
#define DT_IN           prhs[6]

// output arguments
#define MODEL_OUT        plhs[0]
#define SNAPSHOT_OUT    plhs[1]
//#define TEST_OUT        plhs[2] // out argument for test

// the gateway routine
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // begin of declaration
    double *pVelocityModel, *pDataIn, *pModelOut, *pSnapshot;
    double dz, dx, dt;
    int diffOrder, boundary;
    
    int l, i, j, t;
    mwSize nz, nx, nt;
    mwSize pDimsSnapshot[3] = {0};
    
    mxArray *coeff, *oldRtm, *curRtm, *newRtm;
    double *pCoeff, *pOldRtm, *pCurRtm, *pNewRtm;
    
    mxArray *uDampLeft, *vDampLeft, *uDampRight, *vDampRight, *uDampDown, *vDampDown;
    double *puDampLeft, *pvDampLeft, *puDampRight, *pvDampRight, *puDampDown, *pvDampDown;
    
    mxArray *xDampLeft, *xDampRight, *xDamp, *xb, *zDampDown, *zDamp, *zb;
    double *pxDampLeft, *pxDampRight, *pxDamp, *pxb, *pzDampDown, *pzDamp, *pzb;
    
    mxArray *vdtSq;
    double *pVdtSq;
    
    mxArray *source;
    double *pSource;
    
    mxArray *zPhi, *xPhi, *zA, *xA, *zPsi, *xPsi, *zP, *xP;
    double *pzPhi, *pxPhi, *pzA, *pxA, *pzPsi, *pxPsi, *pzP, *pxP;
    
    mxArray *curRtm_diffIn_zPhi, *curRtm_diffOut_zPhi, *curRtm_diffIn_xPhi, *curRtm_diffOut_xPhi;
    double *pCurRtm_diffIn_zPhi, *pCurRtm_diffOut_zPhi, *pCurRtm_diffIn_xPhi, *pCurRtm_diffOut_xPhi;
    
    mxArray *curRtm_diffIn_zA, *curRtm_diffOut_zA, *curRtm_diffIn_xA, *curRtm_diffOut_xA;
    double *pCurRtm_diffIn_zA, *pCurRtm_diffOut_zA, *pCurRtm_diffIn_xA, *pCurRtm_diffOut_xA;
    
    mxArray *zA_diffIn, *zA_diffOut, *xA_diffIn, *xA_diffOut;
    double *pzA_diffIn, *pzA_diffOut, *pxA_diffIn, *pxA_diffOut;
    
    // end of declaration
    
    if (nrhs < 7)
        mexErrMsgTxt("All 7 input arguments shall be provided!");
    
    // ATTENTION: mxGetPr might just produce a 1D array that is linearized according to Matlab convention (column order)
    pVelocityModel = mxGetPr(VM_IN);
    pDataIn = mxGetPr(DATA_IN);
    diffOrder = *mxGetPr(DIFFORDER_IN);
    boundary = *mxGetPr(BOUNDARY_IN);
    dz = *mxGetPr(DZ_IN);
    dx = *mxGetPr(DX_IN);
    dt = *mxGetPr(DT_IN);
    
    nz = mxGetM(VM_IN);
    nx = mxGetN(VM_IN);
    mxAssert(nx == mxGetM(DATA_IN), "Velocity model and input data should have the same x-axis grids!");
    nt = mxGetN(DATA_IN);
    
    /*mexPrintf("pVelocityModel[1] = %f, pDataIn[1] = %f, diffOrder = %d, boundary = %d, dz = %f, dx = %f, dt = %f\nnz = %d, nx = %d, nt = %d\n",
            pVelocityModel[1], pDataIn[1], diffOrder, boundary, dz, dx, dt, nz, nx, nt);*/
    
    // initialize storage    
    pDimsSnapshot[0] = nz;
    pDimsSnapshot[1] = nx;
    pDimsSnapshot[2] = nt;
    SNAPSHOT_OUT = mxCreateNumericArray(3, pDimsSnapshot, mxDOUBLE_CLASS, mxREAL);
    pSnapshot = mxGetPr(SNAPSHOT_OUT);
    
    coeff = dCoef(diffOrder, "s");
    pCoeff = mxGetPr(coeff);
    l = 2 * diffOrder - 1;
    
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
    memcpy(pvDampRight, pVelocityModel + (nx-boundary) * nz, sizeof(double) * nz * boundary);
    xDampRight = dampPml(uDampRight, vDampRight, boundary * dx);
    pxDampRight = mxGetPr(xDampRight);
    
    xDamp = mxCreateDoubleMatrix(nz, nx, mxREAL);
    pxDamp = mxGetPr(xDamp);
    memcpy(pxDamp, pxDampLeft, sizeof(double) * nz * boundary);
    memcpy(pxDamp + (nx-boundary) * nz, pxDampRight, sizeof(double) * nz * boundary);
    
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
    // rtm(:, :, 1) - oldRtm; rtm(:, :, 2) - curRtm; rtm(:, :, 3) - newRtm
    oldRtm = mxCreateDoubleMatrix(nz+2*l, nx+2*l, mxREAL);
    pOldRtm = mxGetPr(oldRtm);
    curRtm = mxCreateDoubleMatrix(nz+2*l, nx+2*l, mxREAL);
    pCurRtm = mxGetPr(curRtm);
    newRtm = mxCreateDoubleMatrix(nz+2*l, nx+2*l, mxREAL);
    pNewRtm = mxGetPr(newRtm);
    
    zPhi = mxCreateDoubleMatrix(nz+2*l, nx, mxREAL);
    pzPhi = mxGetPr(zPhi);
    xPhi = mxCreateDoubleMatrix(nz, nx+2*l, mxREAL);
    pxPhi = mxGetPr(xPhi);
    zA = mxCreateDoubleMatrix(nz+2*l, nx, mxREAL);
    pzA = mxGetPr(zA);
    xA = mxCreateDoubleMatrix(nz, nx+2*l, mxREAL);
    pxA = mxGetPr(xA);
    zPsi = mxCreateDoubleMatrix(nz+l, nx, mxREAL);
    pzPsi = mxGetPr(zPsi);
    xPsi = mxCreateDoubleMatrix(nz, nx+l, mxREAL);
    pxPsi = mxGetPr(xPsi);
    zP = mxCreateDoubleMatrix(nz+l, nx, mxREAL);
    pzP = mxGetPr(zP);
    xP = mxCreateDoubleMatrix(nz, nx+l, mxREAL);
    pxP = mxGetPr(xP);
    
    vdtSq = mxCreateDoubleMatrix(nz, nx, mxREAL);
    pVdtSq = mxGetPr(vdtSq);
    for (j = 0; j < nx; j++)
        for (i = 0; i < nz; i++)
            pVdtSq[j * nz + i] = (pVelocityModel[j * nz + i] * dt) * (pVelocityModel[j * nz + i] * dt);
    
    curRtm_diffIn_zPhi = mxCreateDoubleMatrix(nz+l, nx, mxREAL);
    pCurRtm_diffIn_zPhi = mxGetPr(curRtm_diffIn_zPhi);
    curRtm_diffIn_xPhi = mxCreateDoubleMatrix(nz, nx+l, mxREAL);
    pCurRtm_diffIn_xPhi = mxGetPr(curRtm_diffIn_xPhi);
    curRtm_diffIn_zA = mxCreateDoubleMatrix(nz+2*l, nx, mxREAL);
    pCurRtm_diffIn_zA = mxGetPr(curRtm_diffIn_zA);
    curRtm_diffIn_xA = mxCreateDoubleMatrix(nz, nx+2*l, mxREAL);
    pCurRtm_diffIn_xA = mxGetPr(curRtm_diffIn_xA);
    zA_diffIn = mxCreateDoubleMatrix(nz+l, nx, mxREAL);
    pzA_diffIn = mxGetPr(zA_diffIn);
    xA_diffIn = mxCreateDoubleMatrix(nz, nx+l, mxREAL);
    pxA_diffIn = mxGetPr(xA_diffIn);
    
    // izi = l:(nz+l-1); // len: nz
    // ixi = l:(nx+l-1); // len: nx
    // izl = (diffOrder-1):(nz+2*l-diffOrder-1); // len: nz+l
    // ixl = (diffOrder-1):(nx+2*l-diffOrder-1); // len: nx+l
    for (t = nt-1; t >= 0; t--)     // reverse propagation
    {
        //source = zeros(nz, nx);
        //source(1, :) = data(:, it).';
        source = mxCreateDoubleMatrix(nz, nx, mxREAL);
        pSource = mxGetPr(source);
        for (j = 0; j < nx; j++)
            pSource[j * nz] = pDataIn[t * nx + j];
        
        // zPhi(izi, :) = zb .* zPhi(izi, :) + (zb - 1) .* diffOperator(rtm(izl+1, ixi, 2), coeff, dz, 1);
        for (j = l; j < nx+l; j++)
            for (i = diffOrder; i < nz+2*l-diffOrder+1; i++)
                pCurRtm_diffIn_zPhi[(j - l) * (nz+l) + (i-diffOrder)] = pCurRtm[j * (nz+2*l) + i];
        curRtm_diffOut_zPhi = diffOperator(curRtm_diffIn_zPhi, coeff, dz, 1);
        pCurRtm_diffOut_zPhi = mxGetPr(curRtm_diffOut_zPhi);
        
        for (j = 0; j < nx; j++)
            for (i = l; i < nz + l; i++)
                pzPhi[j * (nz+2*l) + i] = pzb[j * nz + (i - l)] * pzPhi[j * (nz+2*l) + i] +
                        (pzb[j * nz + (i - l)] - 1) * pCurRtm_diffOut_zPhi[j * nz + (i - l)];
        
        // xPhi(:, ixi) = xb .* xPhi(:, ixi) + (xb - 1) .* diffOperator(rtm(izi, ixl+1, 2), coeff, dx, 2);
        for (j = diffOrder; j < nx+2*l-diffOrder+1; j++)
            for (i = l; i < nz+l; i++)
                pCurRtm_diffIn_xPhi[(j-diffOrder) * nz + (i - l)] = pCurRtm[j * (nz+2*l) + i];
        curRtm_diffOut_xPhi = diffOperator(curRtm_diffIn_xPhi, coeff, dx, 2);
        pCurRtm_diffOut_xPhi = mxGetPr(curRtm_diffOut_xPhi);
        
        for (j = l; j < nx + l; j++)
            for (i = 0; i < nz; i++)
                pxPhi[j * nz + i] = pxb[(j - l) * nz + i] * pxPhi[j * nz + i] +
                        (pxb[(j - l) * nz + i] - 1) * pCurRtm_diffOut_xPhi[(j - l) * nz + i];
        
        // zA(izl, :) = diffOperator(rtm(:, ixi, 2), coeff, dz, 1) + zPhi(izl, :);
        memcpy(pCurRtm_diffIn_zA, pCurRtm + l * (nz+2*l), sizeof(double) * nx * (nz+2*l));
        curRtm_diffOut_zA = diffOperator(curRtm_diffIn_zA, coeff, dz, 1);
        pCurRtm_diffOut_zA = mxGetPr(curRtm_diffOut_zA);
        
        for (j = 0; j < nx; j++)
            for (i = diffOrder - 1; i < nz+2*l-diffOrder; i++)
                pzA[j * (nz+2*l) + i] = pCurRtm_diffOut_zA[j * (nz+l) + (i - (diffOrder - 1))] + pzPhi[j * (nz+2*l) + i];
        
        // xA(:, ixl) = diffOperator(rtm(izi, :, 2), coeff, dx, 2) + xPhi(:, ixl);
        for (j = 0; j < nx+2*l; j++)
            for (i = l; i < nz+l; i++)
                pCurRtm_diffIn_xA[j * nz + (i - l)] = pCurRtm[j * (nz+2*l) + i];
        curRtm_diffOut_xA = diffOperator(curRtm_diffIn_xA, coeff, dx, 2);
        pCurRtm_diffOut_xA = mxGetPr(curRtm_diffOut_xA);
        
        for (j = diffOrder - 1; j < nx+2*l-diffOrder; j++)
            for (i = 0; i < nz; i++)
                pxA[j * nz + i] = pCurRtm_diffOut_xA[(j - (diffOrder - 1)) * nz + i] + pxPhi[j * nz + i];
        
        // zPsi(izi, :) = zb .* zPsi(izi, :) + (zb - 1) .* diffOperator(zA(izl, :), coeff, dz, 1);
        for (j = 0; j < nx; j++)
            for (i = diffOrder - 1; i < nz+2*l-diffOrder; i++)
                pzA_diffIn[j * (nz+l) + (i - (diffOrder - 1))] = pzA[j * (nz+2*l) + i];
        zA_diffOut = diffOperator(zA_diffIn, coeff, dz, 1);
        pzA_diffOut = mxGetPr(zA_diffOut);
        
        for (j = 0; j < nx; j++)
            for (i = l; i < nz + l; i++)
                pzPsi[j * (nz+l) + i] = pzb[j * nz + (i - l)] * pzPsi[j * (nz+l) + i] +
                        (pzb[j * nz + (i - l)] - 1) * pzA_diffOut[j * nz + (i - l)];
        
        // xPsi(:, ixi) = xb .* xPsi(:, ixi) + (xb - 1) .* diffOperator(xA(:, ixl), coeff, dx, 2);
        for (j = diffOrder - 1; j < nx+2*l-diffOrder; j++)
            for (i = 0; i < nz; i++)
                pxA_diffIn[(j - (diffOrder - 1)) * nz + i] = pxA[j * nz + i];
        xA_diffOut = diffOperator(xA_diffIn, coeff, dx, 2);
        pxA_diffOut = mxGetPr(xA_diffOut);
        
        for (j = l; j < nx + l; j++)
            for (i = 0; i < nz; i++)
                pxPsi[j * nz + i] = pxb[(j - l) * nz + i] * pxPsi[j * nz + i] +
                        (pxb[(j - l) * nz + i] - 1) * pxA_diffOut[(j - l) * nz + i];
        
        // zP(izi, :) = diffOperator(zA(izl, :), coeff, dz, 1) + zPsi(izi, :);
        for (j = 0; j < nx; j++)
            for (i = l; i < nz + l; i++)
                pzP[j * (nz+l) + i] = pzA_diffOut[j * nz + (i - l)] + pzPsi[j * (nz+l) + i];
        
        // xP(:, ixi) = diffOperator(xA(:, ixl), coeff, dx, 2) + xPsi(:, ixi);
        for (j = l; j < nx + l; j++)
            for (i = 0; i < nz; i++)
                pxP[j * nz + i] = pxA_diffOut[(j - l) * nz + i] + pxPsi[j * nz + i];
        
        /* ======================================================================
         * One-step finite difference calculation
         * ====================================================================== */
        // rtm(izi, ixi, 3) = vdtSq .* (zP(izi, :) + xP(:, ixi) - source) + 2 * rtm(izi, ixi, 2) - rtm(izi, ixi, 1);
        for (j = l; j < nx + l; j++)
            for (i = l; i < nz + l; i++)
                pNewRtm[j * (nz+2*l) + i] = pVdtSq[(j - l) * nz + (i - l)] *
                        ( pzP[(j - l) * (nz+l) + i] + pxP[j * nz + (i - l)] - pSource[(j - l) * nz + (i - l)] ) +
                        2 * pCurRtm[j * (nz+2*l) + i] - pOldRtm[j * (nz+2*l) + i];
        
        // update finite difference matrices
        // rtm(:, :, 1) = rtm(:, :, 2);
        memcpy(pOldRtm, pCurRtm, sizeof(double) * (nz+2*l) * (nx+2*l));
        
        // rtm(:, :, 2) = rtm(:, :, 3);
        memcpy(pCurRtm, pNewRtm, sizeof(double) * (nz+2*l) * (nx+2*l));
        
        // update snapshot
        // snapshot(:, :, it) = rtm(izi, ixi, 2);
        for (j = 0; j < nx; j++)
            for (i = 0; i < nz; i++)
                pSnapshot[t * (nz * nx) + j * nz + i] = pCurRtm[(j + l) * (nz+2*l) + (i + l)];
        
        // ATTENTION: Don't forget to free dynamic memory allocated by MXCREATE* functions (except for output arrays), otherwise memory leak will occur
        mxDestroyArray(source);
        mxDestroyArray(curRtm_diffOut_zPhi);
        mxDestroyArray(curRtm_diffOut_xPhi);
        mxDestroyArray(curRtm_diffOut_zA);
        mxDestroyArray(curRtm_diffOut_xA);
        mxDestroyArray(zA_diffOut);
        mxDestroyArray(xA_diffOut);
    }
    
    // write out final wavefield
    // model = rtm(:, :, 1);
    MODEL_OUT = mxCreateDoubleMatrix(nz+2*l, nx+2*l, mxREAL);
    pModelOut = mxGetPr(MODEL_OUT);
    memcpy(pModelOut, pOldRtm, sizeof(double) * (nz+2*l) * (nx+2*l));
    
    // test begin
    // TEST_OUT = source;
    // test end
    
    // ATTENTION: Don't forget to free dynamic memory allocated by MXCREATE* functions (except for output arrays), otherwise memory leak will occur
    mxDestroyArray(coeff);
    mxDestroyArray(oldRtm);
    mxDestroyArray(curRtm);
    mxDestroyArray(newRtm);
    mxDestroyArray(uDampLeft);
    mxDestroyArray(vDampLeft);
    mxDestroyArray(uDampRight);
    mxDestroyArray(vDampRight);
    mxDestroyArray(uDampDown);
    mxDestroyArray(vDampDown);
    mxDestroyArray(xDampLeft);
    mxDestroyArray(xDampRight);
    mxDestroyArray(xDamp);
    mxDestroyArray(xb);
    mxDestroyArray(zDampDown);
    mxDestroyArray(zDamp);
    mxDestroyArray(zb);
    mxDestroyArray(vdtSq);
    mxDestroyArray(zPhi);
    mxDestroyArray(xPhi);
    mxDestroyArray(zA);
    mxDestroyArray(xA);
    mxDestroyArray(zPsi);
    mxDestroyArray(xPsi);
    mxDestroyArray(zP);
    mxDestroyArray(xP);
    mxDestroyArray(curRtm_diffIn_zPhi);
    mxDestroyArray(curRtm_diffIn_xPhi);
    mxDestroyArray(curRtm_diffIn_zA);
    mxDestroyArray(curRtm_diffIn_xA);
    mxDestroyArray(zA_diffIn);
    mxDestroyArray(xA_diffIn);
}