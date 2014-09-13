/*
  Copyright (C) 2004 Caltech
  Written by Lexing Ying
*/

#include "mex.h"
#include "matrix.h"

#include "fdct_usfft.hpp"

#include "mexaux.hpp"

using namespace std;
using namespace fdct_usfft_ns;

//digital curvelet transform
extern void _main();

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  //mexErrMsgTxt("done");
  if(nrhs!=5)
	 mexErrMsgTxt("5 inputs required");
  if(nlhs!=6)
	 mexErrMsgTxt("6 outputs required");
  
  int N1; mex2cpp(prhs[0], N1);
  int N2; mex2cpp(prhs[1], N2);
  int nbscales; mex2cpp(prhs[2], nbscales);
  int nbangles_coarse; mex2cpp(prhs[3], nbangles_coarse);
  int allcurvelets; mex2cpp(prhs[4], allcurvelets);
  
  vector< vector<dbl2> > sx;  vector< vector<dbl2> > sy;
  vector< vector<double> > fx;  vector< vector<double> > fy;
  vector< vector<int> >    nx;  vector< vector<int> >    ny;
  fdct_usfft_param(N1, N2, nbscales, nbangles_coarse, allcurvelets, sx, sy, fx, fy, nx, ny);
  
  cpp2mex(sx, plhs[0]);
  cpp2mex(sy, plhs[1]);
  cpp2mex(fx, plhs[2]);
  cpp2mex(fy, plhs[3]);
  cpp2mex(nx, plhs[4]);
  cpp2mex(ny, plhs[5]);
    
  return;
}
