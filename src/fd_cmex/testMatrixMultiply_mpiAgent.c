#include "mex.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/* data related */
#define A_IN        prhs[0] /* (M * P) */
#define B_IN        prhs[1] /* (P * N) */
#define C_OUT       plhs[0] /* (M * N) */
#define TASKID_OUT  plhs[1]
/* MPI related */
#define MASTER      0       /* taskid of first task */
#define FROM_MASTER	1       /* self-defined message tag */
#define FROM_WORKER 2       /* self-defined message tag */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* data related declaration */
    double *pA, *pB, *pC;
    int i, j, k;
    mwSize M, P, N;
    
    /* MPI related declaration */
    int errorCode,                  /* MPI error code */
            numTasks,               /* number of tasks in partition */
            taskId,                 /* task identifier */
            numWorkers,             /* number of worker tasks */
            source,                 /* taskid of message source */
            dest,                   /* taskid of message destination */
            msgType,                /* message type */
            cols,                   /* columns of matrix B sent to each worker */
            avgCol, extra, offset;  /* variables used to determine columns of matrix B sent to each worker */
    MPI_Status status;
    
    M = mxGetM(A_IN);               /* rows of matrix A */
    P = mxGetN(A_IN);               /* columns of matrix A and rows of matrix B */
    mxAssert(P == mxGetM(B_IN), "Columns of A must be equal to rows of B!");
    N = mxGetN(B_IN);               /* columns of matrix B */
    
    // test begin
    mexPrintf("output: M = %d, P = %d, N = %d\n", M, P, N);
    // test end
    
    /* initialize the MPI environment */
    errorCode = MPI_Init(NULL, NULL);
    if (errorCode != MPI_SUCCESS)
    {
        mexPrintf("MPI System: Error starting MPI program, terminating...\n");
        MPI_Abort(MPI_COMM_WORLD, errorCode);
    }
    /* find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
    /* find out process rank  */
    MPI_Comm_rank(MPI_COMM_WORLD, &taskId);
    
    if (numTasks < 2)
    {
        mexPrintf("MPI System: Need at least two MPI tasks, terminating...\n");
        MPI_Abort(MPI_COMM_WORLD, errorCode);
    }
    numWorkers = numTasks - 1;
    
    /* master */
    if (taskId == MASTER)
    {
        pA = mxGetPr(A_IN);
        pB = mxGetPr(B_IN);
        pC = (double*)mxCalloc(M * N, sizeof(double));
        // test begin
        mexPrintf("master (test output): pA[0] = %f, pB[0] = %f\n", pA[0], pB[0]);
        // test end
        
        mexPrintf("MPI System: started with %d tasks, including %d workers\n", numTasks, numWorkers);
        avgCol = N / numWorkers;
        extra = N % numWorkers;
        offset = 0;
        /* send data to workers */
        msgType = FROM_MASTER;
        for (dest = 1; dest <= numWorkers; dest++)
        {
            cols = (dest <= extra) ? (avgCol + 1) : (avgCol); /* number of columns processed by current worker */
            MPI_Send(&offset, 1, MPI_INT, dest, msgType, MPI_COMM_WORLD);
            // test begin
            mexPrintf("master (test output): sent offset = %d to worker %d\n", offset, dest);
            // test end
            MPI_Send(&cols, 1, MPI_INT, dest, msgType, MPI_COMM_WORLD);
            // test begin
            mexPrintf("master (test output): sent cols = %d to worker %d\n", cols, dest);
            // test end
            MPI_Send(pA, M * P, MPI_DOUBLE, dest, msgType, MPI_COMM_WORLD);
            // test begin
            mexPrintf("master (test output): sent matrix A (%d * %d) to worker %d\n", M, P, dest);
            mexPrintf("master (test output): pA[0] = %f, pA[%d * %d - 1] = %f sent\n", pA[0], M, P, pA[M * P - 1]);
            // test end
            MPI_Send(pB + P * offset, P * cols, MPI_DOUBLE, dest, msgType, MPI_COMM_WORLD); /* Matlab is using column indexing */
            // test begin
            mexPrintf("master (test output): pB[%d * %d] = %f, pB[%d * %d + %d * %d - 1] = %f sent\n", P, offset, pB[P * offset], P, offset, P, cols, pB[P * offset + P * cols - 1]);
            // test end
            mexPrintf("Master: sent %d columns of matrix B to Worker %d, offset = %d\n", cols, dest, offset);
            offset += cols;
        }
        /* receive results from workers */
        msgType = FROM_WORKER;
        for (source = 1; source <= numWorkers; source++)
        {
            MPI_Recv(&offset, 1, MPI_INT, source, msgType, MPI_COMM_WORLD, &status);
            // test begin
            mexPrintf("master (test output): received offset = %d\n", offset);
            // test end
            MPI_Recv(&cols, 1, MPI_INT, source, msgType, MPI_COMM_WORLD, &status);
            // test begin
            mexPrintf("master (test output): received cols = %d\n", cols);
            // test end
            MPI_Recv(pC + offset * M, M * cols, MPI_DOUBLE, source, msgType, MPI_COMM_WORLD, &status);
            // test begin
            mexPrintf("master (test output): pC[%d * %d] = %f received\n", offset, M, pC[offset * M]);
            // test end
            mexPrintf("Master: received %d columns of result matrix C from Worker %d\n", cols, source);
        }
    }
    /* workers */
    else
    {
        /* receive data from master */
        msgType = FROM_MASTER;
        MPI_Recv(&offset, 1, MPI_INT, MASTER, msgType, MPI_COMM_WORLD, &status);
        // test begin
        mexPrintf("worker %d (test output): received offset = %d from master\n", taskId, offset);
        // test end
        MPI_Recv(&cols, 1, MPI_INT, MASTER, msgType, MPI_COMM_WORLD, &status);
        // test begin
        mexPrintf("worker %d (test output): received cols = %d from master\n", taskId, cols);
        // test end
        pA = (double*)mxCalloc(M * P, sizeof(double));
        MPI_Recv(pA, M * P, MPI_DOUBLE, MASTER, msgType, MPI_COMM_WORLD, &status);
        // test begin
        mexPrintf("worker %d (test output): received matrix A (%d * %d) from master\n", taskId, M, P);
        mexPrintf("worker %d (test output): pA[0] = %f, pA[%d * %d - 1] = %f received\n", taskId, pA[0], M, P, pA[M * P - 1]);
        // test end
        pB = (double*)mxCalloc(P * cols, sizeof(double));
        MPI_Recv(pB, P * cols, MPI_DOUBLE, MASTER, msgType, MPI_COMM_WORLD, &status);
        // test begin
        mexPrintf("worker %d (test output): pB[0] = %f, pB[%d * %d - 1] = %f received\n", taskId, pB[0], P, cols, pB[P * cols - 1]);
        // test end
        mexPrintf("Worker %d: received %d columns of matrix B from master\n", taskId, cols);
        // test begin
        mexPrintf("worker %d (test output): calculating multiplication of (%d * %d) with (%d * %d) matrices\n", taskId, M, P, P, cols);
        // test end
        pC = (double*)mxCalloc(M * cols, sizeof(double));
        for (i = 0; i < M; i++)
            for (j = 0; j < cols; j++)
                for (k = 0; k < P; k++)
                    pC[j * M + i] += pA[k * M + i] * pB[j * P + k];
        // test begin
        mexPrintf("worker %d (test output): pC[0] = %f\n", taskId, pC[0]);
        // test end
        
        /* send results to master */
        msgType = FROM_WORKER;
        MPI_Send(&offset, 1, MPI_INT, MASTER, msgType, MPI_COMM_WORLD);
        MPI_Send(&cols, 1, MPI_INT, MASTER, msgType, MPI_COMM_WORLD);
        MPI_Send(pC, M * cols, MPI_DOUBLE, MASTER, msgType, MPI_COMM_WORLD);
        mexPrintf("Worker %d: sent %d cols of the matrix multiplication result to Master\n", taskId, cols);
        
        /* free memory */
        mxFree(pA);
        mxFree(pB);
        mxFree(pC);
    }
    
    /* shut down MPI */
    MPI_Finalize();
    
    /* output */
    C_OUT = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    if (taskId == MASTER)
    {
        mxSetPr(C_OUT, pC);
        mxSetM(C_OUT, M);
        mxSetN(C_OUT, N);
        // test begin
        mexPrintf("master (taskId = %d) (test output): final size of C_OUT is (%d * %d)\n", taskId, mxGetM(C_OUT), mxGetN(C_OUT));
        // test end
    }
    TASKID_OUT = mxCreateNumericMatrix(1, 1, mxUINT8_CLASS, mxREAL);
    *((int*)mxGetData(TASKID_OUT)) = taskId;
    
    /* exit Matlab */
    /*mexCallMATLAB(0, 0, 0, 0, "exit");*/
}