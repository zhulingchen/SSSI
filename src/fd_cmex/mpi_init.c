/* ======================================================================
 *
 * mpi_init.c
 *
 * a mex wrapper of MPI_Init because it is erroneous to call MPI_Init and
 * MPI_Finalize multiple times
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int flagFinalized, flagInitialized;
    MPI_Finalized(&flagFinalized);
    if (flagFinalized == 0)	/* no MPI has finalized before*/
    {
        MPI_Initialized(&flagInitialized);
        if (flagInitialized == 0)   /* no MPI has initialized before*/
            MPI_Init(NULL, NULL);
        else
        {
            mexPrintf("MPI System: MPI has already initialized!\n");
            MPI_Abort(MPI_COMM_WORLD, flagInitialized);
        }
    }
    else	/* MPI has already finalized before*/
    {
        mexPrintf("MPI System: MPI has already finalized before! It is erroneous to call MPI_Init() after MPI_FINALIZE was invoked.\n");
        MPI_Abort(MPI_COMM_WORLD, flagFinalized);
    }
}