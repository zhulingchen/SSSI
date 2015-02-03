/* ======================================================================
 *
 * fwdTimeCpmlFor2dAw_openmpi_mex.c
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
#define TEST_OUT        plhs[3]     /* out argument for test */

/* MPI-related macros */
#define LEFT_TO_RIGHT   1
#define RIGHT_TO_LEFT   2


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
    double *pTestOut;
    mwSize pDimsTestOut[3] = {0};
    int *sendcounts_band2_nx, *displs_band2_nx;
    /* test end */
    
    /* MPI-related variables */
    int numProcesses, taskId, errorCode;
    int avg_nx, rem_nx, block_nx, offset_block_nx, recvcount_block_nx;
    int block_nx_dampPml;
    int *sendcounts_block_nx, *displs_block_nx, *sendcounts_band_nx, *displs_band_nx;
    MPI_Datatype type_ztPlane, type_ztPlane_resized, type_ztxBlock, type_ztxBlock_resized;
    MPI_Request send_request, recv_request;
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
    // test begin
    sendcounts_band2_nx = (int*)mxCalloc(numProcesses, sizeof(int));
    displs_band2_nx = (int*)mxCalloc(numProcesses, sizeof(int));
    // test end
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
        // test begin
        sendcounts_band2_nx[irank] = (nz+2*l) * block_nx;
        displs_band2_nx[irank] = (nz+2*l) * offset_block_nx;
        // test end
        offset_block_nx += sendcounts_block_nx[irank];
    }
    recvcount_block_nx = sendcounts_block_nx[taskId];
    
    /* scatter velocity model */
    pVelocityModel_local = (double*)mxCalloc(nz * recvcount_block_nx, sizeof(double));
    MPI_Scatterv(pVelocityModel, sendcounts_band_nx, displs_band_nx, MPI_DOUBLE,
            pVelocityModel_local, nz * recvcount_block_nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     // test begin
//     mexPrintf("\nMPI_Scatterv for velocity model Done!\ndispls_block_nx[%d] = %d, sendcounts_block_nx[%d] = %d\n",
//             taskId, displs_block_nx[taskId], taskId, recvcount_block_nx);
//     mexPrintf("worker %d: pVelocityModel_local[%d] = %f, pVelocityModel_local[%d] = %f, pVelocityModel_local[%d] = %f\n",
//             taskId, 0, pVelocityModel_local[0], nz, pVelocityModel_local[nz], nz * recvcount_block_nx - 1, pVelocityModel_local[nz * recvcount_block_nx - 1]);
//     // test end
    
    /* create a strided vector datatype for send */
    if (taskId == 0)
    {
        // send datatype
        MPI_Type_vector(nt, nz, nz*nx, MPI_DOUBLE, &type_ztPlane);
        MPI_Type_commit(&type_ztPlane);
        MPI_Type_create_resized(type_ztPlane, 0, nz * sizeof(double), &type_ztPlane_resized);
        MPI_Type_commit(&type_ztPlane_resized);
    }
    
    /* create a strided vector datatype for receive */
    MPI_Type_vector(nt, nz, nz*recvcount_block_nx, MPI_DOUBLE, &type_ztxBlock);
    MPI_Type_commit(&type_ztxBlock);
    MPI_Type_create_resized(type_ztxBlock, 0, nz * sizeof(double), &type_ztxBlock_resized);
    MPI_Type_commit(&type_ztxBlock_resized);
    
    /* scatter source field */
    pSource_local = (double*)mxCalloc(nz * recvcount_block_nx * nt, sizeof(double));
    MPI_Scatterv(pSource, sendcounts_block_nx, displs_block_nx, type_ztPlane_resized,
            pSource_local, recvcount_block_nx, type_ztxBlock_resized, 0, MPI_COMM_WORLD);
//     // test begin
//     mexPrintf("\nMPI_Scatterv for source field Done!\n");
//     mexPrintf("worker %d: pSource_local[%d] = %f, pSource_local[%d] = %f, pSource_local[%d] = %f, pSource_local[%d] = %f, pSource_local[%d] = %f\n",
//             taskId, 0, pSource_local[0], 1, pSource_local[1],
//             nz * recvcount_block_nx - 1, pSource_local[nz * recvcount_block_nx - 1], nz * recvcount_block_nx, pSource_local[nz * recvcount_block_nx],
//             nz * recvcount_block_nx * nt - 1, pSource_local[nz * recvcount_block_nx * nt - 1]);
//     // test end
    
    /* x-axis damp profile (left), for those tasks whose grids are in left artificial boundary */
    pxDamp_local = (double*)mxCalloc(nz * recvcount_block_nx, sizeof(double));
    if (displs_block_nx[taskId] < boundary)
    {
        if (displs_block_nx[taskId] + recvcount_block_nx <= boundary)
        {
            // test begin
            mexPrintf("worker: %d: all grids in the region are in left artificial boundary\n", taskId);
            // test end
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
            // test begin
            mexPrintf("worker: %d: some grids are in left artificial boundary, and some are not\n", taskId);
            // test end
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
            // test begin
            mexPrintf("worker: %d: some grids are in left artificial boundary, and some are in right artificial boundary\n", taskId);
            // test end
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
            // test begin
            mexPrintf("worker: %d: all grids in the region are in right artificial boundary\n", taskId);
            // test end
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
            // test begin
            mexPrintf("worker: %d: some grids are in right artificial boundary, and some are not\n", taskId);
            // test end
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
    else
    {
        // test begin
        mexPrintf("worker: %d: none grids in the region is in left/right artificial boundary\n", taskId);
        // test end
        /* none grids in the region is in left/right artificial boundary */
    }
    pxb_local = (double*)mxCalloc(nz * recvcount_block_nx, sizeof(double));
    for (j = 0; j < recvcount_block_nx; j++)
        for (i = 0; i < nz; i++)
            pxb_local[j * nz + i] = exp(-pxDamp_local[j * nz + i] * dt);
    mxFree(pxDamp_local);
    
//     // test begin
//     mexPrintf("worker: %d:\npxb_local[%d] = %f, pxb_local[%d] = %f,\npxb_local[%d] = %f, pxb_local[%d] = %f\npxb_local[%d] = %f\n",
//             taskId, 0, pxb_local[0], 1, pxb_local[1],
//             nz-1, pxb_local[nz-1], nz, pxb_local[nz],
//             nz * recvcount_block_nx-1, pxb_local[nz * recvcount_block_nx-1]);
//     // test end
    
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
    
//     // test begin
//     mexPrintf("worker: %d:\npzb_local[%d] = %f, pzb_local[%d] = %f,\npzb_local[%d] = %f, pzb_local[%d] = %f\npzb_local[%d] = %f\n",
//             taskId, 0, pzb_local[0], 1, pzb_local[1],
//             nz-1, pzb_local[nz-1], nz, pzb_local[nz],
//             nz * recvcount_block_nx-1, pzb_local[nz * recvcount_block_nx-1]);
//     mexPrintf("worker: %d:\npVdtSq_local[%d] = %f, pVdtSq_local[%d] = %f,\npVdtSq_local[%d] = %f, pVdtSq_local[%d] = %f\npVdtSq_local[%d] = %f\n",
//             taskId, 0, pVdtSq_local[0], 1, pVdtSq_local[1],
//             nz-1, pVdtSq_local[nz-1], nz, pVdtSq_local[nz],
//             nz * recvcount_block_nx-1, pVdtSq_local[nz * recvcount_block_nx-1]);
//     pTestOut = (double*)mxCalloc(nz * nx, sizeof(double));
//     MPI_Gatherv(pxb_local, nz * recvcount_block_nx, MPI_DOUBLE,
//             pTestOut, sendcounts_band_nx, displs_band_nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     pTestOut = (double*)mxCalloc(nz * nx * nt, sizeof(double));
//     MPI_Gatherv(pSource_local, recvcount_block_nx, type_ztxBlock_resized,
//             pTestOut, sendcounts_block_nx, displs_block_nx, type_ztPlane_resized, 0, MPI_COMM_WORLD);
//     // test end
    
    
    
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
    
    /* finite difference time domain method */
    /*
     * izi = l:(nz+l-1); len: nz
     * ixi = l:(recvcount_block_nx+l-1); len: recvcount_block_nx
     * izl = (diffOrder-1):(nz+2*l-diffOrder-1); len: nz+l
     * ixl = (diffOrder-1):(recvcount_block_nx+2*l-diffOrder-1); len: recvcount_block_nx+l
     */
    for (t = 0; t < nt; t++)
    {
        /* transfer of ghost cell values */
        if ( taskId > 0 )
        {
            MPI_Isend(pCurFdm_local + (nz+2*l) * l, (nz+2*l) * l, MPI_DOUBLE,
                    taskId - 1, RIGHT_TO_LEFT, MPI_COMM_WORLD, &send_request);
        }
        if ( taskId < numProcesses - 1 )
        {
            MPI_Irecv(pCurFdm_local + (nz+2*l) * (nx+l), (nz+2*l) * l, MPI_DOUBLE,
                    taskId + 1, RIGHT_TO_LEFT, MPI_COMM_WORLD, &recv_request);
        }
        if ( taskId < numProcesses - 1 )
        {
            MPI_Isend(pCurFdm_local + (nz+2*l) * nx, (nz+2*l) * l, MPI_DOUBLE,
                    taskId + 1, LEFT_TO_RIGHT, MPI_COMM_WORLD, &send_request);
        }
        if ( taskId > 0 )
        {
            MPI_Irecv(pCurFdm_local, (nz+2*l) * l, MPI_DOUBLE,
                    taskId - 1, LEFT_TO_RIGHT, MPI_COMM_WORLD, &recv_request);
        }
        // test begin
        mexPrintf("NOTICE: t = %d, worker: %d, data communication complete!\n", t, taskId);
        // test end
        /* complete a nonblocking communication */
        MPI_Wait(&send_request, &status);
        MPI_Wait(&recv_request, &status);
        
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
    }
    
    
    // test begin
    pTestOut = (double*)mxCalloc((nz+2*l) * nx, sizeof(double));
    MPI_Gatherv(pzPhi_local, (nz+2*l) * recvcount_block_nx, MPI_DOUBLE,
            pTestOut, sendcounts_band2_nx, displs_band2_nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // test end
    
    
    
    
    /* free datatype */
    if (taskId == 0)
    {
        MPI_Type_free(&type_ztPlane);
        MPI_Type_free(&type_ztPlane_resized);
    }
    MPI_Type_free(&type_ztxBlock);
    MPI_Type_free(&type_ztxBlock_resized);
    /* free dynamic array */
    mxFree(sendcounts_block_nx);
    mxFree(displs_block_nx);
    mxFree(sendcounts_band_nx);
    mxFree(displs_band_nx);
    // test begin
    mxFree(sendcounts_band2_nx);
    mxFree(displs_band2_nx);
    // test end
    mxFree(pVelocityModel_local);
    mxFree(pSource_local);
    
    /* shut down MPI */
    MPI_Finalize();
    
    /* output */
    DATA_OUT = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    SNAPSHOT_OUT = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    /* test begin */
    TEST_OUT = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    /* test end */
    if (taskId == 0)
    {
        /* initialize output storage */
        //DATA_OUT = mxCreateDoubleMatrix(nx, nt, mxREAL);
        //pData = mxGetPr(DATA_OUT);
        pData = (double*)mxCalloc(nx * nt, sizeof(double));
        mxSetPr(DATA_OUT, pData);
        mxSetM(DATA_OUT, nx);
        mxSetN(DATA_OUT, nt);
        
        pSnapshot = (double*)mxCalloc(nz * nx * nt, sizeof(double));
        pDimsSnapshot[0] = nz;
        pDimsSnapshot[1] = nx;
        pDimsSnapshot[2] = nt;
        //SNAPSHOT_OUT = mxCreateNumericArray(3, pDimsSnapshot, mxDOUBLE_CLASS, mxREAL);
        //pSnapshot = mxGetPr(SNAPSHOT_OUT);
        mxSetPr(SNAPSHOT_OUT, pSnapshot);
        mxSetDimensions(SNAPSHOT_OUT, pDimsSnapshot, 3);
        
        /* test begin */
        mxSetPr(TEST_OUT, pTestOut);
        
        mxSetM(TEST_OUT, nz+2*l);
        mxSetN(TEST_OUT, nx);
//         pDimsTestOut[0] = nz;
//         pDimsTestOut[1] = nx;
//         pDimsTestOut[2] = nt;
//         mxSetDimensions(TEST_OUT, pDimsTestOut, 3);
        /* test end */
    }
    TASKID_OUT = mxCreateNumericMatrix(1, 1, mxUINT8_CLASS, mxREAL);
    *((int*)mxGetData(TASKID_OUT)) = taskId;
    
    
    /* free dynamic array */
    mxFree(pCoeff);
}