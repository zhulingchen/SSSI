/* ======================================================================
 *
 * mpi_finalize.c
 *
 * a mex wrapper of MPI_Finalize because it is erroneous to call MPI_Init
 * and MPI_Finalize multiple times
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
    int flagInitialized, flagFinalized;
    MPI_Initialized(&flagInitialized);
    if (flagInitialized == 1)   /* MPI has initialized before */
    {
        MPI_Finalized(&flagFinalized);
        if (flagFinalized == 0) /* no MPI has finalized before*/
            MPI_Finalize();
        else
        {
            mexPrintf("MPI System: MPI has already finalized!\n");
            MPI_Abort(MPI_COMM_WORLD, flagFinalized);
        }
    }
    else    /* no MPI has initialized before */
    {
        mexPrintf("MPI System: No MPI has initialized before and therefore no MPI_Finalize is required!\n");
        MPI_Abort(MPI_COMM_WORLD, flagInitialized);
    }
}