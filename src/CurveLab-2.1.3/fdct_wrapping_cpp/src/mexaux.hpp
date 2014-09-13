/*
   Copyright (C) 2004 Caltech
   Written by Lexing Ying
*/

#ifndef _MEXAUX_HPP_
#define _MEXAUX_HPP_

#include "mex.h"
#include "matrix.h"

#include "fdct_wrapping.hpp"

FDCT_WRAPPING_NS_BEGIN_NAMESPACE

inline void mex2cpp(const mxArray*& md, int& cd);
inline void cpp2mex(const int& cd, mxArray*& md);

inline void mex2cpp(const mxArray*& md, double& cd);
inline void cpp2mex(const double& cd, mxArray*& md);

inline void mex2cpp(const mxArray*& md, CpxOffMat& cd);
inline void cpp2mex(const CpxOffMat& cd, mxArray*& md);

inline void mex2cpp(const mxArray*& md, CpxNumMat& cd);
inline void cpp2mex(const CpxNumMat& cd, mxArray*& md);

template <class T> inline void mex2cpp(const mxArray*& md, vector<T>& cd);
template <class T> inline void cpp2mex(const vector<T>& cd, mxArray*& md);

//----------------------int
inline void mex2cpp(const mxArray*& md, int& cd)
{
  cd = int(mxGetScalar(md));
  return;
}
inline void cpp2mex(const int& cd, mxArray*& md)
{
  md = mxCreateDoubleScalar(cd);
  return;
}

//----------------------double
inline void mex2cpp(const mxArray*& md, double& cd)
{
  cd = mxGetScalar(md);
  return;
}
inline void cpp2mex(const double& cd, mxArray*& md)
{
  md = mxCreateDoubleScalar(cd);
  return;
}

//----------------------cpxoffmat
inline void mex2cpp(const mxArray*& md, CpxOffMat& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md);
  double* xr = mxGetPr(md);
  double* xi = mxGetPi(md);
  int s = -m/2;
  int t = -n/2;
  cd.resize(m,n);
  if(xr!=NULL && xi!=NULL) {
	 int cnt = 0;
	 for(int j=t; j<t+n; j++)
		for(int i=s; i<s+m; i++) {
		  cd(i,j) = cpx(xr[cnt], xi[cnt]);
		  cnt++;
		}
  } else if(xr!=NULL && xi==NULL) {
	 int cnt = 0;
	 for(int j=t; j<t+n; j++)
		for(int i=s; i<s+m; i++) {
		  cd(i,j) = cpx(xr[cnt], 0);
		  cnt++;
		}
  } else if(xr==NULL && xi!=NULL) {
	 int cnt = 0;
	 for(int j=t; j<t+n; j++)
		for(int i=s; i<s+m; i++) {
		  cd(i,j) = cpx(0, xi[cnt]);
		  cnt++;
		}
  }
  return;
}
inline void cpp2mex(const CpxOffMat& cd, mxArray*& md)
{
  int m = cd.m();
  int n = cd.n();
  int s = -m/2;
  int t = -n/2;
  md = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
  double* xr = mxGetPr(md);
  double* xi = mxGetPi(md);
  int cnt = 0;
  for(int j=t; j<t+n; j++)
	 for(int i=s; i<s+m; i++) {
		xr[cnt] = real(cd(i,j));
		xi[cnt] = imag(cd(i,j));
		cnt++;
	 }
  return;
}


//----------------------cpxnummat
inline void mex2cpp(const mxArray*& md, CpxNumMat& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md);
  double* xr = mxGetPr(md);
  double* xi = mxGetPi(md);
  cd.resize(m,n);
  if(xr!=NULL && xi!=NULL) {
	 int cnt = 0;
	 for(int j=0; j<n; j++)
		for(int i=0; i<m; i++) {
		  cd(i,j) = cpx(xr[cnt], xi[cnt]);
		  cnt++;
		}
  } else if(xr!=NULL && xi==NULL) {
	 int cnt = 0;
	 for(int j=0; j<n; j++)
		for(int i=0; i<m; i++) {
		  cd(i,j) = cpx(xr[cnt], 0);
		  cnt++;
		}
  } else if(xr==NULL && xi!=NULL) {
	 int cnt = 0;
	 for(int j=0; j<n; j++)
		for(int i=0; i<m; i++) {
		  cd(i,j) = cpx(0, xi[cnt]);
		  cnt++;
		}
  }
  return;
}
inline void cpp2mex(const CpxNumMat& cd, mxArray*& md)
{
  int m = cd.m();
  int n = cd.n();
  md = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
  double* xr = mxGetPr(md);
  double* xi = mxGetPi(md);
  int cnt = 0;
  for(int j=0; j<n; j++)
	 for(int i=0; i<m; i++) {
		xr[cnt] = real(cd(i,j));
		xi[cnt] = imag(cd(i,j));
		cnt++;
	 }
  return;
}

//----------------------vector<...>
template <class T> inline void mex2cpp(const mxArray*& md, vector<T>& cd)
{
  int m = mxGetM(md); assert(m==1);
  int n = mxGetN(md);
  cd.resize(n);
  for(int ci=0; ci<n; ci++) {
	 const mxArray*tt = mxGetCell(md, ci);
	 mex2cpp(tt, cd[ci]);
  }
  return;
}
template <class T> inline void cpp2mex(const vector<T>& cd, mxArray*& md)
{
  int n = cd.size();
  md = mxCreateCellMatrix(1, n);
  for(int ci=0; ci<n; ci++) {
	 mxArray* ss;	 cpp2mex(cd[ci], ss);
	 mxSetCell(md, ci, ss);
  }
  return;
}

FDCT_WRAPPING_NS_END_NAMESPACE

#endif
