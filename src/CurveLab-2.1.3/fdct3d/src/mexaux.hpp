/* FCT (Fast Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#ifndef _MEXAUX_HPP_
#define _MEXAUX_HPP_

#include "mex.h"
#include "matrix.h"

#include "fdct3d.hpp"

inline void mex2cpp(const mxArray*& md, int& cd);
inline void cpp2mex(const int& cd, mxArray*& md);
inline void mex2cpp(const mxArray*& md, double& cd);
inline void cpp2mex(const double& cd, mxArray*& md);
inline void mex2cpp(const mxArray*& md, CpxOffTns& cd);
inline void cpp2mex(CpxOffTns& cd, mxArray*& md);//inline void cpp2mex(const CpxOffTns& cd, mxArray*& md);
template <class T> inline void mex2cpp(const mxArray*& md, vector<T>& cd);
template <class T> inline void cpp2mex(vector<T>& cd, mxArray*& md);//template <class T> inline void cpp2mex(const vector<T>& cd, mxArray*& md);

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

//----------------------cpxofftns
inline void mex2cpp(const mxArray*& md, CpxOffTns& cd)
{
  const int* dims = mxGetDimensions(md);
  int m = dims[0];
  int n = dims[1];
  int p = dims[2];
  double* xr = mxGetPr(md);
  double* xi = mxGetPi(md);
  int s = -m/2;
  int t = -n/2;
  int u = -p/2;
  cd.resize(m,n,p);
  if(xr!=NULL && xi!=NULL) {
	 int cnt = 0;
	 for(int k=u; k<u+p; k++)
		for(int j=t; j<t+n; j++)
		  for(int i=s; i<s+m; i++) {
			 cd(i,j,k) = cpx(xr[cnt], xi[cnt]);
			 cnt++;
		  }
  } else if(xr!=NULL && xi==NULL) {
	 int cnt = 0;
	 for(int k=u; k<u+p; k++)
		for(int j=t; j<t+n; j++)
		  for(int i=s; i<s+m; i++) {
			 cd(i,j,k) = cpx(xr[cnt], 0);
			 cnt++;
		  }
  } else if(xr==NULL && xi!=NULL) {
	 int cnt = 0;
	 for(int k=u; k<u+p; k++)
		for(int j=t; j<t+n; j++)
		  for(int i=s; i<s+m; i++) {
			 cd(i,j,k) = cpx(0, xi[cnt]);
			 cnt++;
		}
  }
  return;
}
inline void cpp2mex(CpxOffTns& cd, mxArray*& md)
{
  int m = cd.m();
  int n = cd.n();
  int p = cd.p();
  int s = -m/2;
  int t = -n/2;
  int u = -p/2;
  int ndim = 3;
  int dims[3];  dims[0] = m;  dims[1] = n;  dims[2] = p;
  md = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxCOMPLEX);
  double* xr = mxGetPr(md);
  double* xi = mxGetPi(md);
  int cnt = 0;
  for(int k=u; k<u+p; k++)
	 for(int j=t; j<t+n; j++)
		for(int i=s; i<s+m; i++) {
		  xr[cnt] = real(cd(i,j,k));
		  xi[cnt] = imag(cd(i,j,k));
		  cnt++;
		}
  cd.resize(0,0,0);
  return;
}

//----------------------cpxnumtns
inline void mex2cpp(const mxArray*& md, CpxNumTns& cd)
{
  const int* dims = mxGetDimensions(md);
  int m = dims[0];
  int n = dims[1];
  int p = dims[2];
  double* xr = mxGetPr(md);
  double* xi = mxGetPi(md);
  cd.resize(m,n,p);
  if(xr!=NULL && xi!=NULL) {
	 int cnt = 0;
	 for(int k=0; k<p; k++)
		for(int j=0; j<n; j++)
		  for(int i=0; i<m; i++) {
			 cd(i,j,k) = cpx(xr[cnt], xi[cnt]);
			 cnt++;
		  }
  } else if(xr!=NULL && xi==NULL) {
	 int cnt = 0;
	 for(int k=0; k<p; k++)
		for(int j=0; j<n; j++)
		  for(int i=0; i<m; i++) {
			 cd(i,j,k) = cpx(xr[cnt], 0);
			 cnt++;
		  }
  } else if(xr==NULL && xi!=NULL) {
	 int cnt = 0;
	 for(int k=0; k<p; k++)
		for(int j=0; j<n; j++)
		  for(int i=0; i<m; i++) {
			 cd(i,j,k) = cpx(0, xi[cnt]);
			 cnt++;
		  }
  }
  return;
}
inline void cpp2mex(CpxNumTns& cd, mxArray*& md)
{
  int m = cd.m();
  int n = cd.n();
  int p = cd.p();
  int ndim = 3;
  int dims[3];  dims[0] = m;  dims[1] = n;  dims[2] = p;
  md = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxCOMPLEX);
  double* xr = mxGetPr(md);
  double* xi = mxGetPi(md);
  int cnt = 0;
  for(int k=0; k<p; k++)
	 for(int j=0; j<n; j++)
		for(int i=0; i<m; i++) {
		  xr[cnt] = real(cd(i,j,k));
		  xi[cnt] = imag(cd(i,j,k));
		  cnt++;
		}
  cd.resize(0,0,0);
  return;
}

//----------------------vector<...>
template <class T> inline void mex2cpp(const mxArray*& md, vector<T>& cd)
{
  int m = mxGetM(md); 
  int n = mxGetN(md);  assert(n==1);
  cd.resize(m*n);
  for(int ci=0; ci<m*n; ci++) {
	 const mxArray*tt = mxGetCell(md, ci);
	 mex2cpp(tt, cd[ci]);
  }
  return;
}
template <class T> inline void cpp2mex(vector<T>& cd, mxArray*& md)
{
  int n = cd.size();
  md = mxCreateCellMatrix(n, 1);
  for(int ci=0; ci<n; ci++) {
	 mxArray* ss;	 cpp2mex(cd[ci], ss);
	 mxSetCell(md, ci, ss);
  }
  return;
}

#endif
