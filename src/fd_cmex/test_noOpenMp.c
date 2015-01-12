#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#define N       1000000

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    int i, nthreads, tid;
    double a[N], sum_result;
    
    /* Some initializations */
    for (i=0; i < N; i++)
        a[i] = i;
    sum_result = 0.0;
    
    for (i=0; i < N; i++)
    {
        sum_result += a[i];
        //mexPrintf("Added %f, Sum = %f\n", a[i], sum_result);
    }
    
    //mexPrintf("Sum = %f\n", sum_result);
    
}