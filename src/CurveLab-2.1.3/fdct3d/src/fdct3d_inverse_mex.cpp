/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#include "mex.h"
#include "matrix.h"

#include "fdct3d.hpp"

#include "mexaux.hpp"

//inverse digital curvelet transform
extern void _main();

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  if(nrhs!=7)
	 mexErrMsgTxt("7 inputs required");
  if(nlhs!=1)
	 mexErrMsgTxt("1 outputs required");
  
  int m; mex2cpp(prhs[0], m);
  int n; mex2cpp(prhs[1], n);
  int p; mex2cpp(prhs[2], p);
  int nbscales; mex2cpp(prhs[3], nbscales);
  int nbdstz_coarse; mex2cpp(prhs[4], nbdstz_coarse);
  int allcurvelets; mex2cpp(prhs[5], allcurvelets);
  vector< vector<CpxNumTns> > c; mex2cpp(prhs[6], c);
  
  CpxNumTns x;
  fdct3d_inverse(m, n, p, nbscales, nbdstz_coarse, allcurvelets, c, x);
  
  cpp2mex(x, plhs[0]);
  
  return;
}
