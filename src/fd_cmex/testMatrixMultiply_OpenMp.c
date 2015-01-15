#include "mex.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define A_IN    prhs[0] /* (M * P) */
#define B_IN    prhs[1] /* (P * N) */
#define C_OUT   plhs[0] /* (M * N) */

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    double *pA, *pB, *pC;
    int i, j, k;
    mwSize M, P, N;
    
    pA = mxGetPr(A_IN);
    pB = mxGetPr(B_IN);
    M = mxGetM(A_IN);
    P = mxGetN(A_IN);
    mxAssert(P == mxGetM(B_IN), "Columns of A must be equal to rows of B!");
    N = mxGetN(B_IN);
    
    C_OUT = mxCreateDoubleMatrix(M, N, mxREAL);
    pC = mxGetPr(C_OUT);
    
#pragma omp parallel private(i, j, k)
    {
#pragma omp for schedule(dynamic)
        for (i = 0; i < M; i++)
            for (j = 0; j < N; j++)
                for (k = 0; k < P; k++)
                    pC[j * M + i] += pA[k * M + i] * pB[j * P + k];
    }
}