/*
  Copyright (C) 2004 Caltech
  Written by Lexing Ying
*/

#include "mex.h"
#include "matrix.h"

#include "fdct_usfft.hpp"

#include "mexaux.hpp"

//inverse digital curvelet transform
extern void _main();

using namespace std;
using namespace fdct_usfft_ns;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  if(nrhs!=6)
	 mexErrMsgTxt("6 inputs required");
  if(nlhs!=1)
	 mexErrMsgTxt("1 outputs required");

  int m; mex2cpp(prhs[0], m);
  int n; mex2cpp(prhs[1], n);
  int nbscales; mex2cpp(prhs[2], nbscales);
  int nbangles_coarse; mex2cpp(prhs[3], nbangles_coarse);
  int allcurvelets; mex2cpp(prhs[4], allcurvelets);
  vector< vector<CpxNumMat> > c; mex2cpp(prhs[5], c);
  
  CpxNumMat x;  //vector<int> extra;
  ifdct_usfft(m, n, nbscales, nbangles_coarse, allcurvelets, c, x);
  
  cpp2mex(x, plhs[0]);
  
  return;
}

