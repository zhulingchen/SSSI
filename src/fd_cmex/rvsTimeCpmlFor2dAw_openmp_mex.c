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
#include <omp.h>
#include "finiteDifference.h"
#include <math.h>
#include <string.h>

/* input arguments */
#define VM_IN           prhs[0]
#define DATA_IN         prhs[1]
#define DIFFORDER_IN	prhs[2]
#define BOUNDARY_IN     prhs[3]
#define DZ_IN           prhs[4]
#define DX_IN           prhs[5]
#define DT_IN           prhs[6]

/* output arguments */
#define MODEL_OUT       plhs[0]
#define SNAPSHOT_OUT    plhs[1]
/*#define TEST_OUT        plhs[2]*/ /* out argument for test */

/* the gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* begin of declaration */
    double *pVelocityModel, *pData, *pModel, *pSnapshot;
    double dz, dx, dt;
    int diffOrder, boundary;
    
    int l, i, j, t;
    mwSize nz, nx, nt;
    mwSize pDimsSnapshot[3] = {0};
    
    
    double *pCoeff, *pOldRtm, *pCurRtm, *pNewRtm;
    double *puDampLeft, *pvDampLeft, *puDampRight, *pvDampRight, *puDampDown, *pvDampDown;
    double *pxDampLeft, *pxDampRight, *pxDamp, *pxb, *pzDampDown, *pzDamp, *pzb;
    double *pVdtSq;
    double *pSource;
    double *pzPhi, *pxPhi, *pzA, *pxA, *pzPsi, *pxPsi, *pzP, *pxP;
    double *pCurRtm_diffIn_zPhi, *pCurRtm_diffOut_zPhi, *pCurRtm_diffIn_xPhi, *pCurRtm_diffOut_xPhi;
    double *pCurRtm_diffIn_zA, *pCurRtm_diffOut_zA, *pCurRtm_diffIn_xA, *pCurRtm_diffOut_xA;
    double *pzA_diffIn, *pzA_diffOut, *pxA_diffIn, *pxA_diffOut;
    
    /* end of declaration */
    
    if (nrhs < 7)
        mexErrMsgTxt("All 7 input arguments shall be provided!");
    
    /* ATTENTION: mxGetPr might just produce a 1D array that is linearized according to Matlab convention (column order) */
    pVelocityModel = mxGetPr(VM_IN);
    pData = mxGetPr(DATA_IN);
    diffOrder = *mxGetPr(DIFFORDER_IN);
    boundary = *mxGetPr(BOUNDARY_IN);
    dz = *mxGetPr(DZ_IN);
    dx = *mxGetPr(DX_IN);
    dt = *mxGetPr(DT_IN);
    
    nz = mxGetM(VM_IN);
    nx = mxGetN(VM_IN);
    mxAssert(nx == mxGetM(DATA_IN), "Velocity model and input data should have the same x-axis grids!");
    nt = mxGetN(DATA_IN);
    
    /* initialize storage */
    pDimsSnapshot[0] = nz;
    pDimsSnapshot[1] = nx;
    pDimsSnapshot[2] = nt;
    SNAPSHOT_OUT = mxCreateNumericArray(3, pDimsSnapshot, mxDOUBLE_CLASS, mxREAL);
    pSnapshot = mxGetPr(SNAPSHOT_OUT);
    
    pCoeff = dCoef(diffOrder, "s");
    l = 2 * diffOrder - 1;
    
    /* damp profile of x-axis */
    puDampLeft = (double*)mxCalloc(nz * boundary, sizeof(double));
    for (j = 0; j < boundary; j++)
        for (i = 0; i < nz; i++)
            puDampLeft[j * nz + i] = (boundary - j) * dx;
    pvDampLeft = (double*)mxCalloc(nz * boundary, sizeof(double));
    memcpy(pvDampLeft, pVelocityModel, sizeof(double) * nz * boundary);
    pxDampLeft = dampPml(puDampLeft, pvDampLeft, nz, boundary, boundary * dx);
    
    puDampRight = (double*)mxCalloc(nz * boundary, sizeof(double));
    for (j = 0; j < boundary; j++)
        for (i = 0; i < nz; i++)
            puDampRight[j * nz + i] = (j + 1) * dx;
    pvDampRight = (double*)mxCalloc(nz * boundary, sizeof(double));
    memcpy(pvDampRight, pVelocityModel + (nx-boundary) * nz, sizeof(double) * nz * boundary);
    pxDampRight = dampPml(puDampRight, pvDampRight, nz, boundary, boundary * dx);
    
    pxDamp = (double*)mxCalloc(nz * nx, sizeof(double));
    memcpy(pxDamp, pxDampLeft, sizeof(double) * nz * boundary);
    memcpy(pxDamp + (nx-boundary) * nz, pxDampRight, sizeof(double) * nz * boundary);
    
    pxb = (double*)mxCalloc(nz * nx, sizeof(double));
    for (j = 0; j < nx; j++)
        for (i = 0; i < nz; i++)
            pxb[j * nz + i] = exp(-pxDamp[j * nz + i] * dt);
    
    /* damp profile of z-axis */
    puDampDown = (double*)mxCalloc(boundary * nx, sizeof(double));
    for (j = 0; j < nx; j++)
        for(i = 0; i < boundary; i++)
            puDampDown[j * boundary + i] = (i + 1) * dz;
    pvDampDown = (double*)mxCalloc(boundary * nx, sizeof(double));
    for (j = 0; j < nx; j++)
        for(i = 0; i < boundary; i++)
            pvDampDown[j * boundary + i] = pVelocityModel[j * nz + (nz - boundary + i)];
    pzDampDown = dampPml(puDampDown, pvDampDown, boundary, nx, boundary * dz);
    
    pzDamp = (double*)mxCalloc(nz * nx, sizeof(double));
    for (j = 0; j < nx; j++)
        for (i = nz-boundary; i < nz; i++)
            pzDamp[j * nz + i] = pzDampDown[j * boundary + i-(nz-boundary)];
    
    pzb = (double*)mxCalloc(nz * nx, sizeof(double));
    for (j = 0; j < nx; j++)
        for (i = 0; i < nz; i++)
            pzb[j * nz + i] = exp(-pzDamp[j * nz + i] * dt);
    
    /* ======================================================================
     * 2-D Acoustic Wave Forward-Time Modeling
     * ====================================================================== */
    /* additional arrays for storage intermediate results */
    /* rtm(:, :, 1) - oldRtm; rtm(:, :, 2) - curRtm; rtm(:, :, 3) - newRtm */
    pOldRtm = (double*)mxCalloc((nz+2*l) * (nx+2*l), sizeof(double));
    pCurRtm = (double*)mxCalloc((nz+2*l) * (nx+2*l), sizeof(double));
    pNewRtm = (double*)mxCalloc((nz+2*l) * (nx+2*l), sizeof(double));
    
    pzPhi = (double*)mxCalloc((nz+2*l) * nx, sizeof(double));
    pxPhi = (double*)mxCalloc(nz * (nx+2*l), sizeof(double));
    pzA = (double*)mxCalloc((nz+2*l) * nx, sizeof(double));
    pxA = (double*)mxCalloc(nz * (nx+2*l), sizeof(double));
    pzPsi = (double*)mxCalloc((nz+l) * nx, sizeof(double));
    pxPsi = (double*)mxCalloc(nz * (nx+l), sizeof(double));
    pzP = (double*)mxCalloc((nz+l) * nx, sizeof(double));
    pxP = (double*)mxCalloc(nz * (nx+l), sizeof(double));
    
    pVdtSq = (double*)mxCalloc(nz * nx, sizeof(double));
    for (j = 0; j < nx; j++)
        for (i = 0; i < nz; i++)
            pVdtSq[j * nz + i] = (pVelocityModel[j * nz + i] * dt) * (pVelocityModel[j * nz + i] * dt);
    
    pCurRtm_diffIn_zPhi = (double*)mxCalloc((nz+l) * nx, sizeof(double));
    pCurRtm_diffIn_xPhi = (double*)mxCalloc(nz * (nx+l), sizeof(double));
    pCurRtm_diffIn_zA = (double*)mxCalloc((nz+2*l) * nx, sizeof(double));
    pCurRtm_diffIn_xA = (double*)mxCalloc(nz * (nx+2*l), sizeof(double));
    pzA_diffIn = (double*)mxCalloc((nz+l) * nx, sizeof(double));
    pxA_diffIn = (double*)mxCalloc(nz * (nx+l), sizeof(double));
    
    /*
     * izi = l:(nz+l-1); len: nz
     * ixi = l:(nx+l-1); len: nx
     * izl = (diffOrder-1):(nz+2*l-diffOrder-1); len: nz+l
     * ixl = (diffOrder-1):(nx+2*l-diffOrder-1); len: nx+l
     */
    for (t = nt-1; t >= 0; t--)     /* reverse propagation */
    {
        /*source = zeros(nz, nx);*/
        /*source(1, :) = data(:, it).';*/
        pSource = (double*)mxCalloc(nz * nx, sizeof(double));
        for (j = 0; j < nx; j++)
            pSource[j * nz] = pData[t * nx + j];
        
        /* zPhi(izi, :) = zb .* zPhi(izi, :) + (zb - 1) .* diffOperator(rtm(izl+1, ixi, 2), coeff, dz, 1); */
        for (j = l; j < nx+l; j++)
            for (i = diffOrder; i < nz+2*l-diffOrder+1; i++)
                pCurRtm_diffIn_zPhi[(j - l) * (nz+l) + (i-diffOrder)] = pCurRtm[j * (nz+2*l) + i];
        pCurRtm_diffOut_zPhi = diffOperator2d(pCurRtm_diffIn_zPhi, nz+l, nx, pCoeff, diffOrder, dz, 1);
        
        for (j = 0; j < nx; j++)
            for (i = l; i < nz + l; i++)
                pzPhi[j * (nz+2*l) + i] = pzb[j * nz + (i - l)] * pzPhi[j * (nz+2*l) + i] +
                        (pzb[j * nz + (i - l)] - 1) * pCurRtm_diffOut_zPhi[j * nz + (i - l)];
        
        /* xPhi(:, ixi) = xb .* xPhi(:, ixi) + (xb - 1) .* diffOperator(rtm(izi, ixl+1, 2), coeff, dx, 2); */
        for (j = diffOrder; j < nx+2*l-diffOrder+1; j++)
            for (i = l; i < nz+l; i++)
                pCurRtm_diffIn_xPhi[(j-diffOrder) * nz + (i - l)] = pCurRtm[j * (nz+2*l) + i];
        pCurRtm_diffOut_xPhi = diffOperator2d(pCurRtm_diffIn_xPhi, nz, nx+l, pCoeff, diffOrder, dx, 2);
        
        for (j = l; j < nx + l; j++)
            for (i = 0; i < nz; i++)
                pxPhi[j * nz + i] = pxb[(j - l) * nz + i] * pxPhi[j * nz + i] +
                        (pxb[(j - l) * nz + i] - 1) * pCurRtm_diffOut_xPhi[(j - l) * nz + i];
        
        /* zA(izl, :) = diffOperator(rtm(:, ixi, 2), coeff, dz, 1) + zPhi(izl, :); */
        memcpy(pCurRtm_diffIn_zA, pCurRtm + l * (nz+2*l), sizeof(double) * nx * (nz+2*l));
        pCurRtm_diffOut_zA = diffOperator2d(pCurRtm_diffIn_zA, nz+2*l, nx, pCoeff, diffOrder, dz, 1);
        
        for (j = 0; j < nx; j++)
            for (i = diffOrder - 1; i < nz+2*l-diffOrder; i++)
                pzA[j * (nz+2*l) + i] = pCurRtm_diffOut_zA[j * (nz+l) + (i - (diffOrder - 1))] + pzPhi[j * (nz+2*l) + i];
        
        /* xA(:, ixl) = diffOperator(rtm(izi, :, 2), coeff, dx, 2) + xPhi(:, ixl); */
        for (j = 0; j < nx+2*l; j++)
            for (i = l; i < nz+l; i++)
                pCurRtm_diffIn_xA[j * nz + (i - l)] = pCurRtm[j * (nz+2*l) + i];
        pCurRtm_diffOut_xA = diffOperator2d(pCurRtm_diffIn_xA, nz, nx+2*l, pCoeff, diffOrder, dx, 2);
        
        for (j = diffOrder - 1; j < nx+2*l-diffOrder; j++)
            for (i = 0; i < nz; i++)
                pxA[j * nz + i] = pCurRtm_diffOut_xA[(j - (diffOrder - 1)) * nz + i] + pxPhi[j * nz + i];
        
        /* zPsi(izi, :) = zb .* zPsi(izi, :) + (zb - 1) .* diffOperator(zA(izl, :), coeff, dz, 1); */
        for (j = 0; j < nx; j++)
            for (i = diffOrder - 1; i < nz+2*l-diffOrder; i++)
                pzA_diffIn[j * (nz+l) + (i - (diffOrder - 1))] = pzA[j * (nz+2*l) + i];
        pzA_diffOut = diffOperator2d(pzA_diffIn, nz+l, nx, pCoeff, diffOrder, dz, 1);
        
        for (j = 0; j < nx; j++)
            for (i = l; i < nz + l; i++)
                pzPsi[j * (nz+l) + i] = pzb[j * nz + (i - l)] * pzPsi[j * (nz+l) + i] +
                        (pzb[j * nz + (i - l)] - 1) * pzA_diffOut[j * nz + (i - l)];
        
        /* xPsi(:, ixi) = xb .* xPsi(:, ixi) + (xb - 1) .* diffOperator(xA(:, ixl), coeff, dx, 2); */
        memcpy(pxA_diffIn, pxA + (diffOrder - 1) * nz, sizeof(double) * (nx+l) * nz);
        pxA_diffOut = diffOperator2d(pxA_diffIn, nz, nx+l, pCoeff, diffOrder, dx, 2);
        
        for (j = l; j < nx + l; j++)
            for (i = 0; i < nz; i++)
                pxPsi[j * nz + i] = pxb[(j - l) * nz + i] * pxPsi[j * nz + i] +
                        (pxb[(j - l) * nz + i] - 1) * pxA_diffOut[(j - l) * nz + i];
        
        /* zP(izi, :) = diffOperator(zA(izl, :), coeff, dz, 1) + zPsi(izi, :); */
        for (j = 0; j < nx; j++)
            for (i = l; i < nz + l; i++)
                pzP[j * (nz+l) + i] = pzA_diffOut[j * nz + (i - l)] + pzPsi[j * (nz+l) + i];
        
        /* xP(:, ixi) = diffOperator(xA(:, ixl), coeff, dx, 2) + xPsi(:, ixi); */
        for (j = l; j < nx + l; j++)
            for (i = 0; i < nz; i++)
                pxP[j * nz + i] = pxA_diffOut[(j - l) * nz + i] + pxPsi[j * nz + i];
        
        /* ======================================================================
         * One-step finite difference calculation
         * ====================================================================== */
        /* rtm(izi, ixi, 3) = vdtSq .* (zP(izi, :) + xP(:, ixi) + source) + 2 * rtm(izi, ixi, 2) - rtm(izi, ixi, 1); */
#pragma omp parallel private(j, i)
        {
#pragma omp for schedule(static, 8)
                    for (j = l; j < nx + l; j++)
                        for (i = l; i < nz + l; i++)
                            pNewRtm[j * (nz+2*l) + i] = pVdtSq[(j - l) * nz + (i - l)] *
                                    ( pzP[(j - l) * (nz+l) + i] + pxP[j * nz + (i - l)] + pSource[(j - l) * nz + (i - l)] ) +
                                    2 * pCurRtm[j * (nz+2*l) + i] - pOldRtm[j * (nz+2*l) + i];
        }
                
                /* update finite difference matrices */
                /* rtm(:, :, 1) = rtm(:, :, 2); */
                memcpy(pOldRtm, pCurRtm, sizeof(double) * (nz+2*l) * (nx+2*l));
                
                /* rtm(:, :, 2) = rtm(:, :, 3); */
                memcpy(pCurRtm, pNewRtm, sizeof(double) * (nz+2*l) * (nx+2*l));
                
                /* update snapshot */
                /* snapshot(:, :, it) = rtm(izi, ixi, 2); */
                for (j = 0; j < nx; j++)
                    for (i = 0; i < nz; i++)
                        pSnapshot[t * (nz * nx) + j * nz + i] = pCurRtm[(j + l) * (nz+2*l) + (i + l)];
                
                /* ATTENTION: Don't forget to free dynamic memory allocated by mxCalloc function (except for output arrays), otherwise memory leak will occur */
                mxFree(pSource);
                mxFree(pCurRtm_diffOut_zPhi);
                mxFree(pCurRtm_diffOut_xPhi);
                mxFree(pCurRtm_diffOut_zA);
                mxFree(pCurRtm_diffOut_xA);
                mxFree(pzA_diffOut);
                mxFree(pxA_diffOut);
    }
    
    /* write out final wavefield */
    /* model = rtm(:, :, 1); */
    MODEL_OUT = mxCreateDoubleMatrix(nz+2*l, nx+2*l, mxREAL);
    pModel = mxGetPr(MODEL_OUT);
    memcpy(pModel, pOldRtm, sizeof(double) * (nz+2*l) * (nx+2*l));
    
    /* test begin */
    /* TEST_OUT = source; */
    /* test end */
    
    /* ATTENTION: Don't forget to free dynamic memory allocated by mxCalloc function (except for output arrays), otherwise memory leak will occur */
    mxFree(pCoeff);
    mxFree(pOldRtm);
    mxFree(pCurRtm);
    mxFree(pNewRtm);
    mxFree(puDampLeft);
    mxFree(pvDampLeft);
    mxFree(puDampRight);
    mxFree(pvDampRight);
    mxFree(puDampDown);
    mxFree(pvDampDown);
    mxFree(pxDampLeft);
    mxFree(pxDampRight);
    mxFree(pxDamp);
    mxFree(pxb);
    mxFree(pzDampDown);
    mxFree(pzDamp);
    mxFree(pzb);
    mxFree(pVdtSq);
    mxFree(pzPhi);
    mxFree(pxPhi);
    mxFree(pzA);
    mxFree(pxA);
    mxFree(pzPsi);
    mxFree(pxPsi);
    mxFree(pzP);
    mxFree(pxP);
    mxFree(pCurRtm_diffIn_zPhi);
    mxFree(pCurRtm_diffIn_xPhi);
    mxFree(pCurRtm_diffIn_zA);
    mxFree(pCurRtm_diffIn_xA);
    mxFree(pzA_diffIn);
    mxFree(pxA_diffIn);
}