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

/* MPI-related macros */
#define MASTER          0
#define TAG_FROM_MASTER 1
#define TAG_FROM_WORKER 2

/* the gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* begin of declaration */
    /* data variables */
    double *pVelocityModel, *pSource, *pData, *pSnapshot;
    double dz, dx, dt;
    int diffOrder, boundary;
    
    int l, i, j, t;
    mwSize nz, nz_vm, nx, nx_vm, nt;
    const mwSize *pDimsSource;
    mwSize pDimsSnapshot[3] = {0};
    
    double *pCoeff, *pOldFdm, *pCurFdm, *pNewFdm;
    double *puDampLeft, *pvDampLeft, *puDampRight, *pvDampRight, *puDampDown, *pvDampDown;
    double *pxDampLeft, *pxDampRight, *pxDamp, *pxb, *pzDampDown, *pzDamp, *pzb;
    double *pVdtSq;
    double *pzPhi, *pxPhi, *pzA, *pxA, *pzPsi, *pxPsi, *pzP, *pxP;
    double *pCurFdm_diffIn_zPhi, *pCurFdm_diffOut_zPhi, *pCurFdm_diffIn_xPhi, *pCurFdm_diffOut_xPhi;
    double *pCurFdm_diffIn_zA, *pCurFdm_diffOut_zA, *pCurFdm_diffIn_xA, *pCurFdm_diffOut_xA;
    double *pzA_diffIn, *pzA_diffOut, *pxA_diffIn, *pxA_diffOut;
    
    /* MPI-related variables */
    int numProcesses, taskId, errorCode;
    int msgType;
    MPI_Status status;
    /* end of declaration */
    
    if (nrhs < 7)
        mexErrMsgTxt("All 7 input arguments shall be provided!");
    
    diffOrder = *mxGetPr(DIFFORDER_IN);
    boundary = *mxGetPr(BOUNDARY_IN);
    dz = *mxGetPr(DZ_IN);
    dx = *mxGetPr(DX_IN);
    dt = *mxGetPr(DT_IN);
    
    pDimsSource = mxGetDimensions(SOURCE_IN);
    nz = pDimsSource[0];
    nz_vm = nz - boundary;
    nx = pDimsSource[1];
    nx_vm = nx - 2 * boundary;
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
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    /* find out process rank  */
    MPI_Comm_rank(MPI_COMM_WORLD, &taskId);
    /* wait for the root processor to do the I/O */
    MPI_Barrier(MPI_COMM_WORLD);
    
    // test begin
    mexPrintf("worker %d: diffOrder = %d, boundary = %d\ndz = %f, dx = %f, dt = %f\nnz = %d, nx = %d, nt = %d\n", taskId, diffOrder, boundary,
            dz, dx, dt, nz, nx, nt);
    for (i = 0; i < diffOrder; i++)
        mexPrintf("pCoeff[%d] = %f\t", i, pCoeff[i]);
    mexPrintf("\n");
    // test end
    
    /* master */
    if (taskId == MASTER)
    {
        /* load velocity model and source field */
        pVelocityModel = mxGetPr(VM_IN);
        pSource = mxGetPr(SOURCE_IN);
        
        /* initialize storage */
        DATA_OUT = mxCreateDoubleMatrix(nx, nt, mxREAL);
        pData = mxGetPr(DATA_OUT);
        
        pDimsSnapshot[0] = nz;
        pDimsSnapshot[1] = nx;
        pDimsSnapshot[2] = nt;
        SNAPSHOT_OUT = mxCreateNumericArray(3, pDimsSnapshot, mxDOUBLE_CLASS, mxREAL);
        pSnapshot = mxGetPr(SNAPSHOT_OUT);
    }
    
    
    /* shut down MPI */
    MPI_Finalize();
    
    
    
    
    TASKID_OUT = mxCreateNumericMatrix(1, 1, mxUINT8_CLASS, mxREAL);
    *((int*)mxGetData(TASKID_OUT)) = taskId;
}