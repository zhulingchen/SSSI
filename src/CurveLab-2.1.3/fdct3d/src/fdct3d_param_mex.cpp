/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#include "mex.h"
#include "matrix.h"

#include "fdct3d.hpp"

#include "mexaux.hpp"

extern void _main();

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  if(nrhs!=6)
	 mexErrMsgTxt("6 inputs required");
  if(nlhs!=6)
	 mexErrMsgTxt("6 outputs required");
  
  int m; mex2cpp(prhs[0], m);
  int n; mex2cpp(prhs[1], n);
  int p; mex2cpp(prhs[2], p);
  int nbscales; mex2cpp(prhs[3], nbscales);
  int nbdstz_coarse; mex2cpp(prhs[4], nbdstz_coarse);
  int allcurvelets; mex2cpp(prhs[5], allcurvelets);
  
  vector< vector<double> > fxs, fys, fzs;
  vector< vector<int   > > nxs, nys, nzs;
  fdct3d_param(m, n, p, nbscales, nbdstz_coarse, allcurvelets, fxs, fys, fzs, nxs, nys, nzs);
  
  cpp2mex(fxs, plhs[0]);
  cpp2mex(fys, plhs[1]);
  cpp2mex(fzs, plhs[2]);
  cpp2mex(nxs, plhs[3]);
  cpp2mex(nys, plhs[4]);
  cpp2mex(nzs, plhs[5]);
  
  return;
}
