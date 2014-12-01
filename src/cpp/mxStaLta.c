/**
 * Short-time-average of long-time-average ratio calculator
 * 
 * SYNTAX:	ratio = mxStaLta(Signal, STW, LTW)
 *
 * INPUT:
 * 		Signal:	Input signal that needs to take STA/LTA ratio
 *		STW:	Short time window 
 * 		LTW:	Long time window
 *
 * OUTPUT:
 *		ratio:	STA/LTA ratio
 *
 * NOTE:
 *		This function only support single channel signal now. Multiple 
 * 		channel support may be added in the future.
 *
 * AUTHOR: 	Lijun Zhu, CeGP, CSIP, Georgia Institute of Technology
 * CONATACT:	gatechzhu@gmail.com
 * DATE:	Nov. 26, 2014
 * 
 * DISCLAIM:	This software is free to use for academic research. No support 
 *		will be provided for its maintainance and I take no 
 *		responsiblity of any damage it may cause by any kind.
**/

/*#include <stdio.h>*/
#include <math.h>
/*#include <gsl/gsl_statistics.h>*/
/*#include <fftw3.h>*/
#include "mex.h"
/*#include "matrix.h"*/

/**
 * MATLAB pass all scalar, vector, and matrix as a 2D array. Only when passing 
 * array more than two dimensions will we encouter mxGetNumberofDimensions() 
 * ~= 2
**/
#define IS_REAL_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_SCALAR(P) (IS_REAL_2D_FULL_DOUBLE(P) && mxGetNumberOfElements(P) == 1)
#define NELEMS(x) (sizeof(x) / sizeof(x[0]))
#define HEAD(N,W) ((N>W-1) ?  N-W+1: 0)


/* Declare function defined later */
double lijun_stats_variance_with_fixed_mean(double*, int, double);

void mexFunction(int nlhs, mxArray *plhs[],		/* Onput variables */
		int nrhs, const mxArray *prhs[])	/* Input variables */
{
	/* Macros for the output and input arguments */
	#define RATIO		plhs[0]
	#define SIGNAL		prhs[0]
	#define STW		prhs[1]
	#define LTW		prhs[2]

	/* Define intermediate variables */
	double *ratio, *signal, mu, tmp;
	int stw, ltw, M, m;

	/* Input check */
	if(nrhs != 3)	/* Check the number of input arguments */
		mexErrMsgTxt("Wrong number of input arguments.");
	else if(nlhs > 1)
		mexErrMsgTxt("Too many output arguments.");

	if(!IS_REAL_2D_FULL_DOUBLE(SIGNAL))	/* Check SIGNAL */
		mexErrMsgTxt("Signal must be a real full double vector");

	if(mxGetN(SIGNAL) > 1)	/* Check the second dimension of input signal */
		mexErrMsgTxt("Signal must be a column vector");

	if(!IS_REAL_SCALAR(STW))	/* Check STA */
		mexErrMsgTxt("STA must be a real scalar");

	if(!IS_REAL_SCALAR(LTW))	/* Check LTA */	
                mexErrMsgTxt("LTA must be a real scalar");
	
	/* Get short-time window and long-time window */
	stw = mxGetScalar(STW);
	ltw = mxGetScalar(LTW);
	

	/* Get signal data and dimension */
	M = mxGetM(SIGNAL);		/* Get vector length */
	signal = mxGetPr(SIGNAL);	/* Get signal pointer */

	/* Check window length comparing with signal length */
	if(stw > ltw)
		mexErrMsgTxt("STW should be smaller than or equal to LTW");

	if(ltw > M)
		mexErrMsgTxt("LTW needs to be smaller than the signal length");

	/* Create output vector */
	RATIO = mxCreateDoubleMatrix(M, 1, mxREAL);
	ratio = mxGetPr(RATIO);

	/* Loop through the sample set to calculate STA/LTA ratio*/
	for(m=0; m<M; m++) {
		/* Test the window MACRO: printf("The window is from %d to %d\n", HEAD(m, sta), m); */
		/*printf("sta is %f\n", lijun_stats_variance_with_fixed_mean(signal+HEAD(m, stw),stw,0));*/
		/*printf("lta is %f\n", lijun_stats_variance_with_fixed_mean(signal+HEAD(m, ltw),ltw,0));*/
		ratio[m] = lijun_stats_variance_with_fixed_mean(signal+HEAD(m, stw),stw,0)/lijun_stats_variance_with_fixed_mean(signal+HEAD(m, ltw), ltw, 0);
	}
}

/* Custom variance function */
double lijun_stats_variance_with_fixed_mean(double* signal, int n, double mu)
{
	int i;
	double var;

	for(i=0; i<n; i++) 
		var += pow(signal[i]-mu, 2);
	
	var /= n;

	return var;
}
