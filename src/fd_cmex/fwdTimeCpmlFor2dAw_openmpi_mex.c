/* ======================================================================
 *
 * fwdTimeCpmlFor2dAw_openmpi_mex.c
 *
 * Simulates 2-d acoustic wave forward propagation using finite difference
 * in time domain with partial differential equation (PDE) using Message
 * Passing Interface (MPI) implementation
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
#include <mpi.h>
#include <math.h>
#include <string.h>
#include "finiteDifference.h"

/* input arguments */
#define VM_IN           prhs[0]
#define SOURCE_IN       prhs[1]
#define DIFFORDER_IN	prhs[2]
#define BOUNDARY_IN     prhs[3]
#define DZ_IN           prhs[4]
#define DX_IN           prhs[5]
#define DT_IN           prhs[6]

/* output arguments */
#define DATA_OUT        plhs[0]
#define SNAPSHOT_OUT    plhs[1]
#define TASKID_OUT      plhs[2]     /* rank (task ID) of the calling process to return */
/* #define TEST_OUT        plhs[3] */     /* out argument for test */

/* MPI-related macros */
#define RIGHT_TO_LEFT   1
#define LEFT_FROM_RIGHT 1
#define LEFT_TO_RIGHT   2
#define RIGHT_FROM_LEFT 2

/* the gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* begin of declaration */
    /* global variables */
    double *pVelocityModel, *pSource, *pData, *pSnapshot;
    double dz, dx, dt;
    int diffOrder, boundary;
    
    int l, irank, i, j, t;
    mwSize nz, nx, nt;
    const mwSize *pDimsSource;
    mwSize pDimsSnapshot[3] = {0};
    
    double *pCoeff;
    
    /* test begin */
    /*
     * double *pTestOut;
     * mwSize pDimsTestOut[3] = {0};
     * int *sendcounts_band2_nx, *displs_band2_nx;
     */
    /* test end */
    
    /* MPI-related variables */
    int numProcesses, taskId, errorCode;
    int avg_nx, rem_nx, block_nx, offset_block_nx, recvcount_block_nx;
    int rem_l_to_left, rem_l_from_right, rem_l_to_right, rem_l_from_left;
    int cur_l_to_left, cur_l_from_right, cur_l_to_right, cur_l_from_left;
    int offset_l_from_right, offset_l_from_left, leftBound, rightBound;
    int block_nx_dampPml;
    int *sendcounts_block_nx, *displs_block_nx, *sendcounts_band_nx, *displs_band_nx;
    MPI_Datatype type_ztPlane_global, type_ztPlane_global_resized, type_ztPlane_local, type_ztPlane_local_resized;
    MPI_Datatype type_trace_local, type_trace_local_resized, type_trace_global, type_trace_global_resized;
    MPI_Request send_request, recv_request;
    int flag;
    MPI_Status status;
    
    /* local variables */
    double *pVelocityModel_local, *pSource_local;
    double *pOldFdm_local, *pCurFdm_local, *pNewFdm_local;
    double *puDampLeft_local, *pvDampLeft_local, *puDampRight_local, *pvDampRight_local, *puDampDown_local, *pvDampDown_local;
    double *pxDampLeft_local, *pxDampRight_local, *pxDamp_local, *pxb_local, *pzDampDown_local, *pzDamp_local, *pzb_local;
    double *pVdtSq_local;
    double *pCurFdm_diffIn_zPhi_local, *pCurFdm_diffIn_xPhi_local, *pCurFdm_diffIn_zA_local,  *pCurFdm_diffIn_xA_local;
    double *pCurFdm_diffOut_zPhi_local, *pCurFdm_diffOut_xPhi_local, *pCurFdm_diffOut_zA_local, *pCurFdm_diffOut_xA_local;
    double *pzA_diffIn_local, *pxA_diffIn_local;
    double *pzA_diffOut_local, *pxA_diffOut_local;
    double *pzPhi_local, *pxPhi_local, *pzA_local, *pxA_local, *pzPsi_local, *pxPsi_local, *pzP_local, *pxP_local;
    double *pData_local, *pSnapshot_local;
    /* end of declaration */
    
    if (nrhs < 7)
        mexErrMsgTxt("All 7 input arguments shall be provided!");
    
    /* load velocity model and source field */
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
    
    pCoeff = dCoef(diffOrder, "s");
    l = 2 * diffOrder - 1;
    
    /* initialize the MPI environment */
    errorCode = MPI_Init(NULL, NULL);
    if (errorCode != MPI_SUCCESS)
    {
        mexPrintf("MPI System: Error starting MPI program, terminating...\n");
        MPI_Abort(MPI_COMM_WORLD, errorCode);
    }
    /* find out number of processes */
    errorCode = MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    /* find out process rank  */
    MPI_Comm_rank(MPI_COMM_WORLD, &taskId);
    
    if (numProcesses > nx)
    {
        mexPrintf("MPI System: Too much MPI tasks, terminating...\n");
        MPI_Abort(MPI_COMM_WORLD, errorCode);
    }
    
    /* set up send counts and displacements for MPI_Scatterv */
    avg_nx = nx / numProcesses;
    rem_nx = nx % numProcesses;
    offset_block_nx = 0;
    sendcounts_block_nx = (int*)mxCalloc(numProcesses, sizeof(int));
    displs_block_nx = (int*)mxCalloc(numProcesses, sizeof(int));
    sendcounts_band_nx = (int*)mxCalloc(numProcesses, sizeof(int));
    displs_band_nx = (int*)mxCalloc(numProcesses, sizeof(int));
    /* test begin */
    /*
     * sendcounts_band2_nx = (int*)mxCalloc(numProcesses, sizeof(int));
     * displs_band2_nx = (int*)mxCalloc(numProcesses, sizeof(int));
     */
    /* test end */
    for (irank = 0; irank < numProcesses; irank++)
    {
        block_nx = (irank < rem_nx) ? (avg_nx + 1) : (avg_nx);
        /* number of ZT-planes processed by each task */
        sendcounts_block_nx[irank] = block_nx;
        /* number of Z-bands processed by each task */
        sendcounts_band_nx[irank] = nz * block_nx;
        /* displacements (relative to sendbuf) from which to take the outgoing elements to each process */
        displs_block_nx[irank] = offset_block_nx;
        displs_band_nx[irank] = nz * offset_block_nx;
        /* test begin */
        /*
         * sendcounts_band2_nx[irank] = (nz+2*l) * block_nx;
         * displs_band2_nx[irank] = (nz+2*l) * offset_block_nx;
         */
        /* test end */
        offset_block_nx += sendcounts_block_nx[irank];
    }
    recvcount_block_nx = sendcounts_block_nx[taskId];
    
    /* scatter velocity model */
    pVelocityModel_local = (double*)mxCalloc(nz * recvcount_block_nx, sizeof(double));
    MPI_Scatterv(pVelocityModel, sendcounts_band_nx, displs_band_nx, MPI_DOUBLE,
            pVelocityModel_local, nz * recvcount_block_nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (taskId == 0)
    {
        /* create a strided vector datatype for send a 2-D plane on global data */
        MPI_Type_vector(nt, nz, nz*nx, MPI_DOUBLE, &type_ztPlane_global);
        MPI_Type_commit(&type_ztPlane_global);
        MPI_Type_create_resized(type_ztPlane_global, 0, nz * sizeof(double), &type_ztPlane_global_resized);
        MPI_Type_commit(&type_ztPlane_global_resized);
        
        /* create a strided vector datatype for receive a data trace on global data */
        MPI_Type_vector(nt, 1, nx, MPI_DOUBLE, &type_trace_global);
        MPI_Type_commit(&type_trace_global);
        MPI_Type_create_resized(type_trace_global, 0, sizeof(double), &type_trace_global_resized);
        MPI_Type_commit(&type_trace_global_resized);
    }
    
    /* create a strided vector datatype for send a 2-D plane on local data */
    MPI_Type_vector(nt, nz, nz*recvcount_block_nx, MPI_DOUBLE, &type_ztPlane_local);
    MPI_Type_commit(&type_ztPlane_local);
    MPI_Type_create_resized(type_ztPlane_local, 0, nz * sizeof(double), &type_ztPlane_local_resized);
    MPI_Type_commit(&type_ztPlane_local_resized);
    
    /* create a strided vector datatype for receive a data trace on local data */
    MPI_Type_vector(nt, 1, recvcount_block_nx, MPI_DOUBLE, &type_trace_local);
    MPI_Type_commit(&type_trace_local);
    MPI_Type_create_resized(type_trace_local, 0, sizeof(double), &type_trace_local_resized);
    MPI_Type_commit(&type_trace_local_resized);
    
    /* scatter source field */
    pSource_local = (double*)mxCalloc(nz * recvcount_block_nx * nt, sizeof(double));
    MPI_Scatterv(pSource, sendcounts_block_nx, displs_block_nx, type_ztPlane_global_resized,
            pSource_local, recvcount_block_nx, type_ztPlane_local_resized, 0, MPI_COMM_WORLD);
    
    /* x-axis damp profile (left), for those tasks whose grids are in left artificial boundary */
    pxDamp_local = (double*)mxCalloc(nz * recvcount_block_nx, sizeof(double));
    if (displs_block_nx[taskId] < boundary)
    {
        if (displs_block_nx[taskId] + recvcount_block_nx <= boundary)
        {
            /* all grids in the region are in left artificial boundary */
            block_nx_dampPml = recvcount_block_nx;
            puDampLeft_local = (double*)mxCalloc(nz * block_nx_dampPml, sizeof(double));
            for (j = 0; j < block_nx_dampPml; j++)
                for (i = 0; i < nz; i++)
                    puDampLeft_local[j * nz + i] = (boundary - displs_block_nx[taskId] - j) * dx;
            pvDampLeft_local = (double*)mxCalloc(nz * block_nx_dampPml, sizeof(double));
            memcpy(pvDampLeft_local, pVelocityModel_local, sizeof(double) * nz * block_nx_dampPml);
            pxDampLeft_local = dampPml(puDampLeft_local, pvDampLeft_local, nz, block_nx_dampPml, boundary * dx);
            
            memcpy(pxDamp_local, pxDampLeft_local, sizeof(double) * nz * block_nx_dampPml);
            mxFree(puDampLeft_local);
            mxFree(pvDampLeft_local);
            mxFree(pxDampLeft_local);
        }
        else
        {
            /* some grids are in left artificial boundary, and some are not */
            block_nx_dampPml = boundary - displs_block_nx[taskId];
            puDampLeft_local = (double*)mxCalloc(nz * block_nx_dampPml, sizeof(double));
            for (j = 0; j < block_nx_dampPml; j++)
                for (i = 0; i < nz; i++)
                    puDampLeft_local[j * nz + i] = (boundary - displs_block_nx[taskId] - j) * dx;
            pvDampLeft_local = (double*)mxCalloc(nz * block_nx_dampPml, sizeof(double));
            memcpy(pvDampLeft_local, pVelocityModel_local, sizeof(double) * nz * block_nx_dampPml);
            pxDampLeft_local = dampPml(puDampLeft_local, pvDampLeft_local, nz, block_nx_dampPml, boundary * dx);
            
            memcpy(pxDamp_local, pxDampLeft_local, sizeof(double) * nz * block_nx_dampPml);
            mxFree(puDampLeft_local);
            mxFree(pvDampLeft_local);
            mxFree(pxDampLeft_local);
        }
        
        if (displs_block_nx[taskId] + recvcount_block_nx > nx-boundary)
        {
            /* some grids are in left artificial boundary, and some are in right artificial boundary */
            block_nx_dampPml = (displs_block_nx[taskId] + recvcount_block_nx) - (nx-boundary);
            puDampRight_local = (double*)mxCalloc(nz * block_nx_dampPml, sizeof(double));
            for (j = 0; j < block_nx_dampPml; j++)
                for (i = 0; i < nz; i++)
                    puDampRight_local[j * nz + i] = (j + 1) * dx;
            pvDampRight_local = (double*)mxCalloc(nz * block_nx_dampPml, sizeof(double));
            memcpy(pvDampRight_local, pVelocityModel_local + nz * (nx-boundary-displs_block_nx[taskId]), sizeof(double) * nz * block_nx_dampPml);
            pxDampRight_local = dampPml(puDampRight_local, pvDampRight_local, nz, block_nx_dampPml, boundary * dx);
            
            memcpy(pxDamp_local + nz * (nx-boundary-displs_block_nx[taskId]), pxDampRight_local, sizeof(double) * nz * block_nx_dampPml);
            mxFree(puDampRight_local);
            mxFree(pvDampRight_local);
            mxFree(pxDampRight_local);
        }
    }
    /* x-axis damp profile (right), for those tasks whose grids are in right artificial boundary */
    else if (displs_block_nx[taskId] + recvcount_block_nx > nx-boundary)
    {
        if (displs_block_nx[taskId] >= nx-boundary)
        {
            /* all grids in the region are in right artificial boundary */
            block_nx_dampPml = recvcount_block_nx;
            puDampRight_local = (double*)mxCalloc(nz * block_nx_dampPml, sizeof(double));
            for (j = 0; j < block_nx_dampPml; j++)
                for (i = 0; i < nz; i++)
                    puDampRight_local[j * nz + i] = (displs_block_nx[taskId] - (nx-boundary) + j + 1) * dx;
            pvDampRight_local = (double*)mxCalloc(nz * block_nx_dampPml, sizeof(double));
            memcpy(pvDampRight_local, pVelocityModel_local, sizeof(double) * nz * block_nx_dampPml);
            pxDampRight_local = dampPml(puDampRight_local, pvDampRight_local, nz, block_nx_dampPml, boundary * dx);
            
            memcpy(pxDamp_local, pxDampRight_local, sizeof(double) * nz * block_nx_dampPml);
            mxFree(puDampRight_local);
            mxFree(pvDampRight_local);
            mxFree(pxDampRight_local);
        }
        else
        {
            /* some grids are in right artificial boundary, and some are not */
            block_nx_dampPml = (displs_block_nx[taskId] + recvcount_block_nx) - (nx-boundary);
            puDampRight_local = (double*)mxCalloc(nz * block_nx_dampPml, sizeof(double));
            for (j = 0; j < block_nx_dampPml; j++)
                for (i = 0; i < nz; i++)
                    puDampRight_local[j * nz + i] = (j + 1) * dx;
            pvDampRight_local = (double*)mxCalloc(nz * block_nx_dampPml, sizeof(double));
            memcpy(pvDampRight_local, pVelocityModel_local + nz * (nx-boundary-displs_block_nx[taskId]), sizeof(double) * nz * block_nx_dampPml);
            pxDampRight_local = dampPml(puDampRight_local, pvDampRight_local, nz, block_nx_dampPml, boundary * dx);
            
            memcpy(pxDamp_local + nz * (nx-boundary-displs_block_nx[taskId]), pxDampRight_local, sizeof(double) * nz * block_nx_dampPml);
            mxFree(puDampRight_local);
            mxFree(pvDampRight_local);
            mxFree(pxDampRight_local);
        }
    }
    /*else*/ /* none grids in the region is in left/right artificial boundary */
    
    pxb_local = (double*)mxCalloc(nz * recvcount_block_nx, sizeof(double));
    for (j = 0; j < recvcount_block_nx; j++)
        for (i = 0; i < nz; i++)
            pxb_local[j * nz + i] = exp(-pxDamp_local[j * nz + i] * dt);
    mxFree(pxDamp_local);
    
    /* z-axis damp profile */
    puDampDown_local = (double*)mxCalloc(boundary * recvcount_block_nx, sizeof(double));
    for (j = 0; j < recvcount_block_nx; j++)
        for(i = 0; i < boundary; i++)
            puDampDown_local[j * boundary + i] = (i + 1) * dz;
    pvDampDown_local = (double*)mxCalloc(boundary * recvcount_block_nx, sizeof(double));
    for (j = 0; j < recvcount_block_nx; j++)
        for(i = 0; i < boundary; i++)
            pvDampDown_local[j * boundary + i] = pVelocityModel_local[j * nz + (nz - boundary + i)];
    pzDampDown_local = dampPml(puDampDown_local, pvDampDown_local, boundary, recvcount_block_nx, boundary * dz);
    
    pzDamp_local = (double*)mxCalloc(nz * recvcount_block_nx, sizeof(double));
    for (j = 0; j < recvcount_block_nx; j++)
        for (i = nz-boundary; i < nz; i++)
            pzDamp_local[j * nz + i] = pzDampDown_local[j * boundary + i-(nz-boundary)];
    
    pzb_local = (double*)mxCalloc(nz * recvcount_block_nx, sizeof(double));
    for (j = 0; j < recvcount_block_nx; j++)
        for (i = 0; i < nz; i++)
            pzb_local[j * nz + i] = exp(-pzDamp_local[j * nz + i] * dt);
    
    mxFree(puDampDown_local);
    mxFree(pvDampDown_local);
    mxFree(pzDampDown_local);
    mxFree(pzDamp_local);
    
    pVdtSq_local = (double*)mxCalloc(nz * recvcount_block_nx, sizeof(double));
    for (j = 0; j < recvcount_block_nx; j++)
        for (i = 0; i < nz; i++)
            pVdtSq_local[j * nz + i] = (pVelocityModel_local[j * nz + i] * dt) * (pVelocityModel_local[j * nz + i] * dt);
    
    /* test begin */
    /*
     * pTestOut = (double*)mxCalloc(nz * nx, sizeof(double));
     * MPI_Gatherv(pxb_local, nz * recvcount_block_nx, MPI_DOUBLE,
     * pTestOut, sendcounts_band_nx, displs_band_nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     * pTestOut = (double*)mxCalloc(nz * nx * nt, sizeof(double));
     * MPI_Gatherv(pSource_local, recvcount_block_nx, type_ztPlane_local_resized,
     * pTestOut, sendcounts_block_nx, displs_block_nx, type_ztPlane_global_resized, 0, MPI_COMM_WORLD);
     */
    /* test end */
    
    /* ======================================================================
     * 2-D Acoustic Wave Forward-Time Modeling
     * ====================================================================== */
    /* additional arrays for storage intermediate results */
    /* oldFdm - old finite difference matrix
     * curFdm - current finite difference matrix
     * newFdm - new finite difference matrix
     */
    pOldFdm_local = (double*)mxCalloc((nz+2*l) * (recvcount_block_nx+2*l), sizeof(double));
    pCurFdm_local = (double*)mxCalloc((nz+2*l) * (recvcount_block_nx+2*l), sizeof(double));
    pNewFdm_local = (double*)mxCalloc((nz+2*l) * (recvcount_block_nx+2*l), sizeof(double));
    
    pzPhi_local = (double*)mxCalloc((nz+2*l) * recvcount_block_nx, sizeof(double));
    pxPhi_local = (double*)mxCalloc(nz * (recvcount_block_nx+2*l), sizeof(double));
    pzA_local = (double*)mxCalloc((nz+2*l) * recvcount_block_nx, sizeof(double));
    pxA_local = (double*)mxCalloc(nz * (recvcount_block_nx+2*l), sizeof(double));
    pzPsi_local = (double*)mxCalloc((nz+l) * recvcount_block_nx, sizeof(double));
    pxPsi_local = (double*)mxCalloc(nz * (recvcount_block_nx+l), sizeof(double));
    pzP_local = (double*)mxCalloc((nz+l) * recvcount_block_nx, sizeof(double));
    pxP_local = (double*)mxCalloc(nz * (recvcount_block_nx+l), sizeof(double));
    
    pCurFdm_diffIn_zPhi_local = (double*)mxCalloc((nz+l) * recvcount_block_nx, sizeof(double));
    pCurFdm_diffIn_xPhi_local = (double*)mxCalloc(nz * (recvcount_block_nx+l), sizeof(double));
    pCurFdm_diffIn_zA_local = (double*)mxCalloc((nz+2*l) * recvcount_block_nx, sizeof(double));
    pCurFdm_diffIn_xA_local = (double*)mxCalloc(nz * (recvcount_block_nx+2*l), sizeof(double));
    pzA_diffIn_local = (double*)mxCalloc((nz+l) * recvcount_block_nx, sizeof(double));
    pxA_diffIn_local = (double*)mxCalloc(nz * (recvcount_block_nx+l), sizeof(double));
    
    pData_local = (double*)mxCalloc(recvcount_block_nx * nt, sizeof(double));
    pSnapshot_local = (double*)mxCalloc(nz * recvcount_block_nx * nt, sizeof(double));
    
    /* finite difference time domain method */
    /*
     * izi = l:(nz+l-1); len: nz
     * ixi = l:(recvcount_block_nx+l-1); len: recvcount_block_nx
     * izl = (diffOrder-1):(nz+2*l-diffOrder-1); len: nz+l
     * ixl = (diffOrder-1):(recvcount_block_nx+2*l-diffOrder-1); len: recvcount_block_nx+l
     */
    for (t = 0; t < nt; t++)
    {
        /* finite difference calculation with PML ABC */
        /* zPhi(izi, :) = zb .* zPhi(izi, :) + (zb - 1) .* diffOperator(fdm(izl+1, ixi, 2), coeff, dz, 1); */
        for (j = l; j < recvcount_block_nx+l; j++)
            for (i = diffOrder; i < nz+2*l-diffOrder+1; i++)
                pCurFdm_diffIn_zPhi_local[(j - l) * (nz+l) + (i-diffOrder)] = pCurFdm_local[j * (nz+2*l) + i];
        pCurFdm_diffOut_zPhi_local = diffOperator2d(pCurFdm_diffIn_zPhi_local, nz+l, recvcount_block_nx, pCoeff, diffOrder, dz, 1);
        
        for (j = 0; j < recvcount_block_nx; j++)
            for (i = l; i < nz + l; i++)
                pzPhi_local[j * (nz+2*l) + i] = pzb_local[j * nz + (i - l)] * pzPhi_local[j * (nz+2*l) + i] +
                        (pzb_local[j * nz + (i - l)] - 1) * pCurFdm_diffOut_zPhi_local[j * nz + (i - l)];
        
        mxFree(pCurFdm_diffOut_zPhi_local);
        
        /* xPhi(:, ixi) = xb .* xPhi(:, ixi) + (xb - 1) .* diffOperator(fdm(izi, ixl+1, 2), coeff, dx, 2); */
        for (j = diffOrder; j < recvcount_block_nx+2*l-diffOrder+1; j++)
            for (i = l; i < nz+l; i++)
                pCurFdm_diffIn_xPhi_local[(j-diffOrder) * nz + (i - l)] = pCurFdm_local[j * (nz+2*l) + i];
        pCurFdm_diffOut_xPhi_local = diffOperator2d(pCurFdm_diffIn_xPhi_local, nz, recvcount_block_nx+l, pCoeff, diffOrder, dx, 2);
        
        for (j = l; j < recvcount_block_nx + l; j++)
            for (i = 0; i < nz; i++)
                pxPhi_local[j * nz + i] = pxb_local[(j - l) * nz + i] * pxPhi_local[j * nz + i] +
                        (pxb_local[(j - l) * nz + i] - 1) * pCurFdm_diffOut_xPhi_local[(j - l) * nz + i];
        
        mxFree(pCurFdm_diffOut_xPhi_local);
        
        /* transfer of ghost cell values */
        if (numProcesses > 1)
        {
            if ((taskId == rem_nx) && (recvcount_block_nx < l))
                rem_l_to_left = l-1;
            else
                rem_l_to_left = l;
            rem_l_from_right = l;
            if ((taskId == rem_nx - 1) && (recvcount_block_nx < l))
                rem_l_to_right = l+1;
            else
                rem_l_to_right = l;
            rem_l_from_left = l;
            offset_l_from_right = 0;
            offset_l_from_left = 0;
            leftBound = 0;
            rightBound = numProcesses - 1;
            
            while ( (rem_l_to_left > 0) || (rem_l_from_right > 0) || (rem_l_to_right > 0) || (rem_l_from_left > 0) )
            {
                mexPrintf("t = %d (pxPhi_local): taskId: %d, rem_l_to_left = %d, rem_l_from_right = %d, rem_l_to_right = %d, rem_l_from_left = %d, leftBound = %d\n", t, taskId, rem_l_to_left, rem_l_from_right, rem_l_to_right, rem_l_from_left, leftBound);
                if ( (taskId > leftBound) && (rem_l_to_left > 0) )
                {
                    cur_l_to_left = (rem_l_to_left <= recvcount_block_nx) ? rem_l_to_left : recvcount_block_nx;
                    MPI_Isend(&cur_l_to_left, 1, MPI_INT, taskId - (leftBound + 1), RIGHT_TO_LEFT, MPI_COMM_WORLD, &send_request);
                    MPI_Wait(&send_request, &status);
                    /* test begin */
                    MPI_Test(&send_request, &flag, &status);
                    mexPrintf("t = %d (pxPhi_local): %d sends to %d, cur_l_to_left(%d) send_request flag: %d\n", t, taskId, taskId - (leftBound + 1), cur_l_to_left, flag);
                    /* test end */
                    MPI_Isend(pxPhi_local + nz * l, nz * cur_l_to_left, MPI_DOUBLE,
                            taskId - (leftBound + 1), RIGHT_TO_LEFT, MPI_COMM_WORLD, &send_request);
                    MPI_Wait(&send_request, &status);
                    /* test begin */
                    MPI_Test(&send_request, &flag, &status);
                    mexPrintf("t = %d (pxPhi_local): %d sends to %d, pxPhi_local send_request flag: %d\n", t, taskId, taskId - (leftBound + 1), flag);
                    /* test end */
                    rem_l_to_left -= cur_l_to_left;
                }
                else
                {
                    rem_l_to_left = 0;
                }
                
                if ( (taskId < rightBound) && (rem_l_from_right > 0) )
                {
                    MPI_Irecv(&cur_l_from_right, 1, MPI_INT, taskId + (leftBound + 1), LEFT_FROM_RIGHT, MPI_COMM_WORLD, &recv_request);
                    MPI_Wait(&recv_request, &status);
                    /* test begin */
                    MPI_Test(&recv_request, &flag, &status);
                    mexPrintf("t = %d (pxPhi_local): %d receives from %d, cur_l_from_right(%d) recv_request flag: %d\n", t, taskId, taskId + (leftBound + 1), cur_l_from_right, flag);
                    /* test end */
                    MPI_Irecv(pxPhi_local + nz * (recvcount_block_nx+l + offset_l_from_right), nz * cur_l_from_right, MPI_DOUBLE,
                            taskId + (leftBound + 1), LEFT_FROM_RIGHT, MPI_COMM_WORLD, &recv_request);
                    MPI_Wait(&recv_request, &status);
                    /* test begin */
                    MPI_Test(&recv_request, &flag, &status);
                    mexPrintf("t = %d (pxPhi_local): %d receives from %d, pxPhi_local recv_request flag: %d\n", t, taskId, taskId + (leftBound + 1), flag);
                    /* test end */
                    rem_l_from_right -= cur_l_from_right;
                    offset_l_from_right += cur_l_from_right;
                }
                else
                {
                    rem_l_from_right = 0;
                }
                
                if ( (taskId < rightBound) && (rem_l_to_right > 0) )
                {
                    cur_l_to_right = (rem_l_to_right <= recvcount_block_nx) ? rem_l_to_right : recvcount_block_nx;
                    MPI_Isend(&cur_l_to_right, 1, MPI_INT, taskId + (leftBound + 1), LEFT_TO_RIGHT, MPI_COMM_WORLD, &send_request);
                    MPI_Wait(&send_request, &status);
                    /* test begin */
                    MPI_Test(&send_request, &flag, &status);
                    mexPrintf("t = %d (pxPhi_local): %d sends to %d, cur_l_to_right(%d) send_request flag: %d\n", t, taskId, taskId + (leftBound + 1), cur_l_to_right, flag);
                    /* test end */
                    MPI_Isend(pxPhi_local + nz * (l + recvcount_block_nx - cur_l_to_right), nz * cur_l_to_right, MPI_DOUBLE,
                            taskId + (leftBound + 1), LEFT_TO_RIGHT, MPI_COMM_WORLD, &send_request);
                    MPI_Wait(&send_request, &status);
                    /* test begin */
                    MPI_Test(&send_request, &flag, &status);
                    mexPrintf("t = %d (pxPhi_local): %d sends to %d, pxPhi_local send_request flag: %d\n", t, taskId, taskId + (leftBound + 1), flag);
                    /* test end */
                    rem_l_to_right -= cur_l_to_right;
                }
                else
                {
                    rem_l_to_right = 0;
                }
                
                if ( (taskId > leftBound) && (rem_l_from_left > 0) )
                {
                    MPI_Irecv(&cur_l_from_left, 1, MPI_INT, taskId - (leftBound + 1), RIGHT_FROM_LEFT, MPI_COMM_WORLD, &recv_request);
                    MPI_Wait(&recv_request, &status);
                    /* test begin */
                    MPI_Test(&recv_request, &flag, &status);
                    mexPrintf("t = %d (pxPhi_local): %d receives from %d, cur_l_from_left(%d) recv_request flag: %d\n", t, taskId, taskId - (leftBound + 1), cur_l_from_left, flag);
                    /* test end */
                    MPI_Irecv(pxPhi_local + nz * (l - cur_l_from_left - offset_l_from_left), nz * cur_l_from_left, MPI_DOUBLE,
                            taskId - (leftBound + 1), RIGHT_FROM_LEFT, MPI_COMM_WORLD, &recv_request);
                    MPI_Wait(&recv_request, &status);
                    /* test begin */
                    MPI_Test(&recv_request, &flag, &status);
                    mexPrintf("t = %d (pxPhi_local): %d receives from %d, pxPhi_local recv_request flag: %d\n", t, taskId, taskId - (leftBound + 1), flag);
                    /* test end */
                    rem_l_from_left -= cur_l_from_left;
                    offset_l_from_left += cur_l_from_left;
                }
                else
                {
                    rem_l_from_left = 0;
                }
                
                leftBound++;
                rightBound--;
            }
        }
        
        mexPrintf("!!!!!!!!!!(pxPhi_local Finished) t = %d, taskId = %d, rem_l_to_left = %d, rem_l_from_right = %d, rem_l_to_right = %d, rem_l_from_left = %d, leftBound = %d, rightBound = %d\n",
                t, taskId, rem_l_to_left, rem_l_from_right, rem_l_to_right, rem_l_from_left, leftBound, rightBound);
        
        /* zA(izl, :) = diffOperator(fdm(:, ixi, 2), coeff, dz, 1) + zPhi(izl, :); */
        memcpy(pCurFdm_diffIn_zA_local, pCurFdm_local + l * (nz+2*l), sizeof(double) * recvcount_block_nx * (nz+2*l));
        pCurFdm_diffOut_zA_local = diffOperator2d(pCurFdm_diffIn_zA_local, nz+2*l, recvcount_block_nx, pCoeff, diffOrder, dz, 1);
        
        for (j = 0; j < recvcount_block_nx; j++)
            for (i = diffOrder - 1; i < nz+2*l-diffOrder; i++)
                pzA_local[j * (nz+2*l) + i] = pCurFdm_diffOut_zA_local[j * (nz+l) + (i - (diffOrder - 1))] + pzPhi_local[j * (nz+2*l) + i];
        
        mxFree(pCurFdm_diffOut_zA_local);
        
        /* xA(:, ixl) = diffOperator(fdm(izi, :, 2), coeff, dx, 2) + xPhi(:, ixl); */
        for (j = 0; j < recvcount_block_nx+2*l; j++)
            for (i = l; i < nz+l; i++)
                pCurFdm_diffIn_xA_local[j * nz + (i - l)] = pCurFdm_local[j * (nz+2*l) + i];
        pCurFdm_diffOut_xA_local = diffOperator2d(pCurFdm_diffIn_xA_local, nz, recvcount_block_nx+2*l, pCoeff, diffOrder, dx, 2);
        
        for (j = diffOrder - 1; j < recvcount_block_nx+2*l-diffOrder; j++)
            for (i = 0; i < nz; i++)
                pxA_local[j * nz + i] = pCurFdm_diffOut_xA_local[(j - (diffOrder - 1)) * nz + i] + pxPhi_local[j * nz + i];
        
        mxFree(pCurFdm_diffOut_xA_local);
        
        /* zPsi(izi, :) = zb .* zPsi(izi, :) + (zb - 1) .* diffOperator(zA(izl, :), coeff, dz, 1); */
        for (j = 0; j < recvcount_block_nx; j++)
            for (i = diffOrder - 1; i < nz+2*l-diffOrder; i++)
                pzA_diffIn_local[j * (nz+l) + (i - (diffOrder - 1))] = pzA_local[j * (nz+2*l) + i];
        pzA_diffOut_local = diffOperator2d(pzA_diffIn_local, nz+l, recvcount_block_nx, pCoeff, diffOrder, dz, 1);
        
        for (j = 0; j < recvcount_block_nx; j++)
            for (i = l; i < nz + l; i++)
                pzPsi_local[j * (nz+l) + i] = pzb_local[j * nz + (i - l)] * pzPsi_local[j * (nz+l) + i] +
                        (pzb_local[j * nz + (i - l)] - 1) * pzA_diffOut_local[j * nz + (i - l)];
        
        /* xPsi(:, ixi) = xb .* xPsi(:, ixi) + (xb - 1) .* diffOperator(xA(:, ixl), coeff, dx, 2); */
        memcpy(pxA_diffIn_local, pxA_local + (diffOrder - 1) * nz, sizeof(double) * (recvcount_block_nx+l) * nz);
        pxA_diffOut_local = diffOperator2d(pxA_diffIn_local, nz, recvcount_block_nx+l, pCoeff, diffOrder, dx, 2);
        
        for (j = l; j < recvcount_block_nx + l; j++)
            for (i = 0; i < nz; i++)
                pxPsi_local[j * nz + i] = pxb_local[(j - l) * nz + i] * pxPsi_local[j * nz + i] +
                        (pxb_local[(j - l) * nz + i] - 1) * pxA_diffOut_local[(j - l) * nz + i];
        
        /* zP(izi, :) = diffOperator(zA(izl, :), coeff, dz, 1) + zPsi(izi, :); */
        for (j = 0; j < recvcount_block_nx; j++)
            for (i = l; i < nz + l; i++)
                pzP_local[j * (nz+l) + i] = pzA_diffOut_local[j * nz + (i - l)] + pzPsi_local[j * (nz+l) + i];
        
        mxFree(pzA_diffOut_local);
        
        /* xP(:, ixi) = diffOperator(xA(:, ixl), coeff, dx, 2) + xPsi(:, ixi); */
        for (j = l; j < recvcount_block_nx + l; j++)
            for (i = 0; i < nz; i++)
                pxP_local[j * nz + i] = pxA_diffOut_local[(j - l) * nz + i] + pxPsi_local[j * nz + i];
        
        mxFree(pxA_diffOut_local);
        
        /* ======================================================================
         * One-step finite difference calculation
         * ====================================================================== */
        /* fdm(izi, ixi, 3) = vdtSq .* (zP(izi, :) + xP(:, ixi) - source(:, :, it)) + 2 * fdm(izi, ixi, 2) - fdm(izi, ixi, 1); */
        for (j = l; j < recvcount_block_nx + l; j++)
            for (i = l; i < nz + l; i++)
                pNewFdm_local[j * (nz+2*l) + i] = pVdtSq_local[(j - l) * nz + (i - l)] *
                        ( pzP_local[(j - l) * (nz+l) + i] + pxP_local[j * nz + (i - l)] - pSource_local[t * (nz * recvcount_block_nx) + (j - l) * nz + (i - l)] ) +
                        2 * pCurFdm_local[j * (nz+2*l) + i] - pOldFdm_local[j * (nz+2*l) + i];
        
        /* update finite difference matrices */
        /* fdm(:, :, 1) = fdm(:, :, 2); */
        memcpy(pOldFdm_local, pCurFdm_local, sizeof(double) * (nz+2*l) * (recvcount_block_nx+2*l));
        
        /* fdm(:, :, 2) = fdm(:, :, 3); */
        memcpy(pCurFdm_local, pNewFdm_local, sizeof(double) * (nz+2*l) * (recvcount_block_nx+2*l));
        
        /* update data */
        /* data(:, it) = fdm(l, ixi, 2); */
        for (i = 0; i < recvcount_block_nx; i++)
            pData_local[t * recvcount_block_nx + i] = pCurFdm_local[(i + l) * (nz+2*l) + l];
        
        /* update snapshot */
        /* snapshot(:, :, it) = fdm(izi, ixi, 2); */
        for (j = 0; j < recvcount_block_nx; j++)
            for (i = 0; i < nz; i++)
                pSnapshot_local[t * (nz * recvcount_block_nx) + j * nz + i] = pCurFdm_local[(j + l) * (nz+2*l) + (i + l)];
        
        /* transfer of ghost cell values */
        if (numProcesses > 1)
        {
            if ((taskId == rem_nx) && (recvcount_block_nx < l))
                rem_l_to_left = l-1;
            else
                rem_l_to_left = l;
            rem_l_from_right = l;
            if ((taskId == rem_nx - 1) && (recvcount_block_nx < l))
                rem_l_to_right = l+1;
            else
                rem_l_to_right = l;
            rem_l_from_left = l;
            offset_l_from_right = 0;
            offset_l_from_left = 0;
            leftBound = 0;
            rightBound = numProcesses - 1;
            
            while ( (rem_l_to_left > 0) || (rem_l_from_right > 0) || (rem_l_to_right > 0) || (rem_l_from_left > 0) )
            {
                mexPrintf("t = %d (pCurFdm_local): taskId: %d, rem_l_to_left = %d, rem_l_from_right = %d, rem_l_to_right = %d, rem_l_from_left = %d, leftBound = %d\n", t, taskId, rem_l_to_left, rem_l_from_right, rem_l_to_right, rem_l_from_left, leftBound);
                if ( (taskId > leftBound) && (rem_l_to_left > 0) )
                {
                    cur_l_to_left = (rem_l_to_left <= recvcount_block_nx) ? rem_l_to_left : recvcount_block_nx;
                    MPI_Isend(&cur_l_to_left, 1, MPI_INT, taskId - (leftBound + 1), RIGHT_TO_LEFT, MPI_COMM_WORLD, &send_request);
                    MPI_Wait(&send_request, &status);
                    /* test begin */
                    MPI_Test(&send_request, &flag, &status);
                    mexPrintf("t = %d (pCurFdm_local): %d sends to %d, cur_l_to_left(%d) send_request flag: %d\n", t, taskId, taskId - (leftBound + 1), cur_l_to_left, flag);
                    /* test end */
                    MPI_Isend(pCurFdm_local + (nz+2*l) * l, (nz+2*l) * cur_l_to_left, MPI_DOUBLE,
                            taskId - (leftBound + 1), RIGHT_TO_LEFT, MPI_COMM_WORLD, &send_request);
                    MPI_Wait(&send_request, &status);
                    /* test begin */
                    MPI_Test(&send_request, &flag, &status);
                    mexPrintf("t = %d (pCurFdm_local): %d sends to %d, pCurFdm_local send_request flag: %d\n", t, taskId, taskId - (leftBound + 1), flag);
                    /* test end */
                    rem_l_to_left -= cur_l_to_left;
                }
                else
                {
                    rem_l_to_left = 0;
                }
                
                if ( (taskId < rightBound) && (rem_l_from_right > 0) )
                {
                    MPI_Irecv(&cur_l_from_right, 1, MPI_INT, taskId + (leftBound + 1), LEFT_FROM_RIGHT, MPI_COMM_WORLD, &recv_request);
                    MPI_Wait(&recv_request, &status);
                    /* test begin */
                    MPI_Test(&recv_request, &flag, &status);
                    mexPrintf("t = %d (pCurFdm_local): %d receives from %d, cur_l_from_right(%d) recv_request flag: %d\n", t, taskId, taskId + (leftBound + 1), cur_l_from_right, flag);
                    /* test end */
                    MPI_Irecv(pCurFdm_local + (nz+2*l) * (recvcount_block_nx+l + offset_l_from_right), (nz+2*l) * cur_l_from_right, MPI_DOUBLE,
                            taskId + (leftBound + 1), LEFT_FROM_RIGHT, MPI_COMM_WORLD, &recv_request);
                    MPI_Wait(&recv_request, &status);
                    /* test begin */
                    MPI_Test(&recv_request, &flag, &status);
                    mexPrintf("t = %d (pCurFdm_local): %d receives from %d, pCurFdm_local recv_request flag: %d\n", t, taskId, taskId + (leftBound + 1), flag);
                    /* test end */
                    rem_l_from_right -= cur_l_from_right;
                    offset_l_from_right += cur_l_from_right;
                }
                else
                {
                    rem_l_from_right = 0;
                }
                
                if ( (taskId < rightBound) && (rem_l_to_right > 0) )
                {
                    cur_l_to_right = (rem_l_to_right <= recvcount_block_nx) ? rem_l_to_right : recvcount_block_nx;
                    MPI_Isend(&cur_l_to_right, 1, MPI_INT, taskId + (leftBound + 1), LEFT_TO_RIGHT, MPI_COMM_WORLD, &send_request);
                    MPI_Wait(&send_request, &status);
                    /* test begin */
                    MPI_Test(&send_request, &flag, &status);
                    mexPrintf("t = %d (pCurFdm_local): %d sends to %d, cur_l_to_right(%d) send_request flag: %d\n", t, taskId, taskId + (leftBound + 1), cur_l_to_right, flag);
                    /* test end */
                    MPI_Isend(pCurFdm_local + (nz+2*l) * (l + recvcount_block_nx - cur_l_to_right), (nz+2*l) * cur_l_to_right, MPI_DOUBLE,
                            taskId + (leftBound + 1), LEFT_TO_RIGHT, MPI_COMM_WORLD, &send_request);
                    MPI_Wait(&send_request, &status);
                    /* test begin */
                    MPI_Test(&send_request, &flag, &status);
                    mexPrintf("t = %d (pCurFdm_local): %d sends to %d, pCurFdm_local send_request flag: %d\n", t, taskId, taskId + (leftBound + 1), flag);
                    /* test end */
                    rem_l_to_right -= cur_l_to_right;
                }
                else
                {
                    rem_l_to_right = 0;
                }
                
                if ( (taskId > leftBound) && (rem_l_from_left > 0) )
                {
                    MPI_Irecv(&cur_l_from_left, 1, MPI_INT, taskId - (leftBound + 1), RIGHT_FROM_LEFT, MPI_COMM_WORLD, &recv_request);
                    MPI_Wait(&recv_request, &status);
//                     /* test begin */
//                     MPI_Test(&recv_request, &flag, &status);
//                     mexPrintf("t = %d (pCurFdm_local): %d receives from %d, cur_l_from_left(%d) recv_request flag: %d\n", t, taskId, taskId - (leftBound + 1), cur_l_from_left, flag);
//                     /* test end */
                    MPI_Irecv(pCurFdm_local + (nz+2*l) * (l - cur_l_from_left - offset_l_from_left), (nz+2*l) * cur_l_from_left, MPI_DOUBLE,
                            taskId - (leftBound + 1), RIGHT_FROM_LEFT, MPI_COMM_WORLD, &recv_request);
                    MPI_Wait(&recv_request, &status);
//                     /* test begin */
//                     MPI_Test(&recv_request, &flag, &status);
//                     mexPrintf("t = %d (pCurFdm_local): %d receives from %d, pCurFdm_local recv_request flag: %d\n", t, taskId, taskId - (leftBound + 1), flag);
//                     /* test end */
                    rem_l_from_left -= cur_l_from_left;
                    offset_l_from_left += cur_l_from_left;
                }
                else
                {
                    rem_l_from_left = 0;
                }
                
                leftBound++;
                rightBound--;
            }
        }
        
        mexPrintf("!!!!!!!!!!(pCurFdm_local Finished) t = %d, taskId = %d, rem_l_to_left = %d, rem_l_from_right = %d, rem_l_to_right = %d, rem_l_from_left = %d, leftBound = %d, rightBound = %d\n",
                t, taskId, rem_l_to_left, rem_l_from_right, rem_l_to_right, rem_l_from_left, leftBound, rightBound);
        
    }
    
    mexPrintf("===== Finished Time Iteration ===== taskId: %d\n", taskId);
    
    /* test begin */
    /*
     * pTestOut = (double*)mxCalloc((nz+2*l) * nx, sizeof(double));
     * MPI_Gatherv(pzPhi_local, (nz+2*l) * recvcount_block_nx, MPI_DOUBLE,
     * pTestOut, sendcounts_band2_nx, displs_band2_nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     */
    /* test end */
    
    /* gather output */
    if (taskId == 0)
    {
        pData = (double*)mxCalloc(nx * nt, sizeof(double));
        pSnapshot = (double*)mxCalloc(nz * nx * nt, sizeof(double));
    }
    MPI_Gatherv(pData_local, recvcount_block_nx, type_trace_local_resized,
            pData, sendcounts_block_nx, displs_block_nx, type_trace_global_resized, 0, MPI_COMM_WORLD);
    MPI_Gatherv(pSnapshot_local, recvcount_block_nx, type_ztPlane_local_resized,
            pSnapshot, sendcounts_block_nx, displs_block_nx, type_ztPlane_global_resized, 0, MPI_COMM_WORLD);
    
    mexPrintf("===== Finished Data Gathering ===== taskId: %d\n", taskId);
    
    /* free datatype */
    if (taskId == 0)
    {
        MPI_Type_free(&type_ztPlane_global);
        MPI_Type_free(&type_ztPlane_global_resized);
        MPI_Type_free(&type_trace_global);
        MPI_Type_free(&type_trace_global_resized);
    }
    MPI_Type_free(&type_ztPlane_local);
    MPI_Type_free(&type_ztPlane_local_resized);
    MPI_Type_free(&type_trace_local);
    MPI_Type_free(&type_trace_local_resized);
    /* free dynamic array */
    mxFree(sendcounts_block_nx);
    mxFree(displs_block_nx);
    mxFree(sendcounts_band_nx);
    mxFree(displs_band_nx);
    mxFree(pVelocityModel_local);
    mxFree(pSource_local);
    mxFree(pOldFdm_local);
    mxFree(pCurFdm_local);
    mxFree(pNewFdm_local);
    mxFree(pxb_local);
    mxFree(pzb_local);
    mxFree(pVdtSq_local);
    mxFree(pzPhi_local);
    mxFree(pxPhi_local);
    mxFree(pzA_local);
    mxFree(pxA_local);
    mxFree(pzPsi_local);
    mxFree(pxPsi_local);
    mxFree(pzP_local);
    mxFree(pxP_local);
    mxFree(pCurFdm_diffIn_zPhi_local);
    mxFree(pCurFdm_diffIn_xPhi_local);
    mxFree(pCurFdm_diffIn_zA_local);
    mxFree(pCurFdm_diffIn_xA_local);
    mxFree(pzA_diffIn_local);
    mxFree(pxA_diffIn_local);
    mxFree(pData_local);
    mxFree(pSnapshot_local);
    
    /* shut down MPI */
    MPI_Finalize();
    
    /* output arrays */
    DATA_OUT = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    SNAPSHOT_OUT = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    /* test begin */
    /* TEST_OUT = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL); */
    /* test end */
    /* set output arrays and dimensions */
    if (taskId == 0)
    {
        mxSetPr(DATA_OUT, pData);
        mxSetM(DATA_OUT, nx);
        mxSetN(DATA_OUT, nt);
        
        pDimsSnapshot[0] = nz;
        pDimsSnapshot[1] = nx;
        pDimsSnapshot[2] = nt;
        mxSetPr(SNAPSHOT_OUT, pSnapshot);
        mxSetDimensions(SNAPSHOT_OUT, pDimsSnapshot, 3);
        
        
        /* test begin */
        /*
         * mxSetPr(TEST_OUT, pTestOut);
         * mxSetM(TEST_OUT, nz);
         * mxSetN(TEST_OUT, nx);
         * pDimsTestOut[0] = nz;
         * pDimsTestOut[1] = nx;
         * pDimsTestOut[2] = nt;
         * mxSetDimensions(TEST_OUT, pDimsTestOut, 3);
         */
        /* test end */
    }
    TASKID_OUT = mxCreateNumericMatrix(1, 1, mxUINT8_CLASS, mxREAL);
    *((int*)mxGetData(TASKID_OUT)) = taskId;
    
    /* free dynamic array */
    mxFree(pCoeff);
}