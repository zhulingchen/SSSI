/*
  Copyright (C) 2004 Caltech
  Written by Lexing Ying
*/

#ifndef _FDCT_USFFT_INLINE_HPP_
#define _FDCT_USFFT_INLINE_HPP_

#include "fdct_usfft.hpp"

FDCT_USFFT_NS_BEGIN_NAMESPACE

inline int fdct_usfft_fftshift(CpxNumMat& T, CpxOffMat& O)
{
  int N1 = T.m();  int N2 = T.n();  O.resize(N1,N2);
  int F1 = N1/2;  int F2 = N2/2; 
  IntOffVec p1(N1);
  for(int i=-F1; i<-F1+N1; i++)	 p1(i) = (i+N1)%N1;
  IntOffVec p2(N2);
  for(int i=-F2; i<-F2+N2; i++)	 p2(i) = (i+N2)%N2;
  for(int j=-F2; j<-F2+N2; j++)
	 for(int i=-F1; i<-F1+N1; i++)
		O(i,j) = T(p1(i),p2(j));
  return 0;
}

inline int fdct_usfft_ifftshift(CpxOffMat& O, CpxNumMat& T)
{
  int N1 = O.m();  int N2 = O.n();  T.resize(N1,N2);
  int F1 = N1/2;  int F2 = N2/2; 
  IntOffVec p1(N1);
  for(int i=-F1; i<-F1+N1; i++)	 p1(i) = (i+N1)%N1;
  IntOffVec p2(N2);
  for(int i=-F2; i<-F2+N2; i++)	 p2(i) = (i+N2)%N2;
  for(int j=-F2; j<-F2+N2; j++)
	 for(int i=-F1; i<-F1+N1; i++)
		T(p1(i),p2(j)) = O(i,j);
  return 0;
}

inline int fdct_usfft_window(double x, double& l, double& r)
{
  if(x<=0) {
	 l = 0;	 r = 1;
  } else if(x>=1) {
	 l = 1;	 r = 0;
  } else {
	 l = exp(1-1/(1-exp(1-1/(1-x))));
	 r = exp(1-1/(1-exp(1-1/x)));
	 double norm = sqrt(l*l+r*r);
	 l /= norm;
	 r /= norm;
  }
  return 0;
}

inline int fdct_usfft_rangecompute(double XL1, double XL2, int& XS1, int& XS2, int& XF1, int& XF2, double& XR1, double& XR2)
{
  XS1 = 2*int(floor(XL1/2))+1;  XS2 = 2*int(floor(XL2/2))+1; //number of samples
  XF1 = int(floor(XL1/2));  XF2 = int(floor(XL2/2)); //offset on either side
  XR1 = XL1/2;  XR2 = XL2/2;
  return 0;
}

inline int fdct_usfft_lowpasscompute(double XL1, double XL2, DblOffMat& lowpass)
{
  int XS1 = 2*int(floor(XL1/2))+1;  int XS2 = 2*int(floor(XL2/2))+1; //number of samples
  int XF1 = int(floor(XL1/2));  int XF2 = int(floor(XL2/2)); //offset on either side
  double XR1 = XL1/2;  double XR2 = XL2/2;
  
  DblOffVec lowpass1(XS1); setvalue(lowpass1,1.0);
  for(int i=-XF1; i<-XR1/2; i++) {
	 double x = (i+XR1)/(XR1/2);
	 double l,r; fdct_usfft_window(x, l, r);
	 lowpass1(i) = l;
	 lowpass1(-i) = l;
  }
  //cerr<<lowpass1;
  DblOffVec lowpass2(XS2); setvalue(lowpass2,1.0);
  for(int i=-XF2; i<-XR2/2; i++) {
	 double x = (i+XR2)/(XR2/2);
	 double l,r; fdct_usfft_window(x, l, r);
	 lowpass2(i) = l;
	 lowpass2(-i) = l;
  }
  //cerr<<lowpass2;
  for(int i=-XF1; i<-XF1+XS1; i++)
	 for(int j=-XF2; j<-XF2+XS2; j++)
		lowpass(i,j) = lowpass1(i) * lowpass2(j);
  return 0;
}

inline int fdct_usfft_fftshift(CpxNumVec& T, CpxOffVec& O)
{
  int N = T.m();  int F = N/2;
  for(int i=-F; i<-F+N; i++) {
	 int s = (i+N)%N;	 //O(i).re = T(s).re;	 O(i).im = T(s).im;
	 O(i) = T(s);
  }
  return 0;
}

inline int fdct_usfft_ifftshift(CpxOffVec& O, CpxNumVec& T)
{
  int N = O.m();  int F = N/2;
  for(int i=-F; i<-F+N; i++) {
	 int s = (i+N)%N;	 //T(s).re = O(i).re;	 T(s).im = O(i).im;
	 T(s) = O(i);
  }
  return 0;
}

inline int fdct_usfft_rotate_forward(int f, double XL1, double XL2, double& TL1, double& TL2)
{
  if(f==0) { //time 1
	 TL1 = XL1;	 TL2 = XL2;
  } else if(f==1) { //time -i 
	 TL1 = XL2;	 TL2 = -XL1;
  } else if(f==2) { //time -1
	 TL1 = -XL1;	 TL2 = -XL2;
  } else if(f==3) { //time i
	 TL1 = -XL2;	 TL2 = XL1;
  }
  return 0;
}

inline int fdct_usfft_rotate_backward(int f, double XL1, double XL2, double& TL1, double& TL2)
{
  if(f==0) {
	 TL1 = XL1;	 TL2 = XL2;
  } else if(f==1) {
	 TL1 = -XL2;	 TL2 = XL1;
  } else if(f==2) {
	 TL1 = -XL1;	 TL2 = -XL2;
  } else if(f==3) {
	 TL1 = XL2;	 TL2 = -XL1;
  }  
  return 0;
}

template <class F>//inline int fdct_usfft_rotate_forward(int f, CpxOffMat& X, CpxOffMat& T)
inline int fdct_usfft_rotate_forward(int f, OffMat<F>& X, OffMat<F>& T)
{
  //rotate face f to the position of face 0
  int XS1 = X.m();  int XS2 = X.n();
  int XF1 =-X.s();  int XF2 =-X.t();
  int TS1, TS2;
  int TF1, TF2;
  if(f==0) {
	 TS1 = XS1;	 TS2 = XS2;
	 TF1 = XF1;	 TF2 = XF2;	 //TL1 = XL1;	 TL2 = XL2;
	 T.resize(TS1, TS2);
	 for(int j=-TF2; j<-TF2+TS2; j++)
		for(int i=-TF1; i<-TF1+TS1; i++)
		  T(i,j) = X(i,j);
  } else if(f==1) {
	 TS1 = XS2;	 TS2 = XS1;
	 TF1 = XF2;	 TF2 = XF1;	 //TL1 = XL2;	 TL2 = XL1;
	 T.resize(TS1, TS2);
	 for(int j=-TF2; j<-TF2+TS2; j++)
		for(int i=-TF1; i<-TF1+TS1; i++)
		  T(i,j) = X(-j,i);
  } else if(f==2) {
	 TS1 = XS1;	 TS2 = XS2;
	 TF1 = XF1;	 TF2 = XF2;	 //TL1 = XL1;	 TL2 = XL2;
	 T.resize(TS1, TS2);
	 for(int j=-TF2; j<-TF2+TS2; j++)
		for(int i=-TF1; i<-TF1+TS1; i++)
		  T(i,j) = X(-i,-j);
  } else if(f==3) {
	 TS1 = XS2;	 TS2 = XS1;
	 TF1 = XF2;	 TF2 = XF1;	 //TL1 = XL2;	 TL2 = XL1;
	 T.resize(TS1, TS2);
	 for(int j=-TF2; j<-TF2+TS2; j++)
		for(int i=-TF1; i<-TF1+TS1; i++)
		  T(i,j) = X(j,-i);
  }
  return 0;
}

template <class F>//inline int fdct_usfft_rotate_backward(int f, CpxOffMat& X, CpxOffMat& T)
inline int fdct_usfft_rotate_backward(int f, OffMat<F>& X, OffMat<F>& T)
{
  int XS1 = X.m();  int XS2 = X.n();
  int XF1 =-X.s();  int XF2 =-X.t();
  int TS1, TS2;
  int TF1, TF2;
  if(f==0) {
	 TS1 = XS1;	 TS2 = XS2;
	 TF1 = XF1;	 TF2 = XF2;	 //TL1 = XL1;	 TL2 = XL2;
	 T.resize(TS1, TS2);
	 for(int j=-TF2; j<-TF2+TS2; j++)
		for(int i=-TF1; i<-TF1+TS1; i++)
		  T(i,j) = X(i,j);
  } else if(f==1) {
	 TS1 = XS2;	 TS2 = XS1;
	 TF1 = XF2;	 TF2 = XF1;	 //TL1 = XL2;	 TL2 = XL1;
	 T.resize(TS1, TS2);
	 for(int j=-TF2; j<-TF2+TS2; j++)
		for(int i=-TF1; i<-TF1+TS1; i++)
		  T(i,j) = X(j,-i);
  } else if(f==2) {
	 TS1 = XS1;	 TS2 = XS2;
	 TF1 = XF1;	 TF2 = XF2;	 //TL1 = XL1;	 TL2 = XL2;
	 T.resize(TS1, TS2);
	 for(int j=-TF2; j<-TF2+TS2; j++)
		for(int i=-TF1; i<-TF1+TS1; i++)
		  T(i,j) = X(-i,-j);
  } else if(f==3) {
	 TS1 = XS2;	 TS2 = XS1;
	 TF1 = XF2;	 TF2 = XF1;	 //TL1 = XL2;	 TL2 = XL1;
	 T.resize(TS1, TS2);
	 for(int j=-TF2; j<-TF2+TS2; j++)
		for(int i=-TF1; i<-TF1+TS1; i++)
		  T(i,j) = X(-j,i);
  }
  return 0;
}

inline double energy(CpxOffMat& m)
{
  double val=0;
  cpx* data = m.data();
  for(int i=0; i<m.m()*m.n(); i++) {	 //val += data[i].re*data[i].re + data[i].im*data[i].im;
	 val += norm(data[i]);
  }
  return val;
}
inline double energy(CpxNumMat& m)
{
  double val=0;
  cpx* data = m.data();
  for(int i=0; i<m.m()*m.n(); i++) {	 //val += data[i].re*data[i].re + data[i].im*data[i].im;
	 val += norm(data[i]);
  }
  return val;
}
inline double energy(DblOffMat& m)
{
  double val=0;
  double* data = m.data();
  for(int i=0; i<m.m()*m.n(); i++) {
	 val += data[i]*data[i];
  }
  return val;
}
inline double energy(DblNumMat& m)
{
  double val=0;
  double* data = m.data();
  for(int i=0; i<m.m()*m.n(); i++) {
	 val += data[i]*data[i];
  }
  return val;
}

inline double energy(CpxOffVec& m)
{
  double val=0;
  cpx* data = m.data();
  for(int i=0; i<m.m(); i++) {	 //val += data[i].re*data[i].re + data[i].im*data[i].im;
	 val += norm(data[i]);
  }
  return val;
}
inline double energy(CpxNumVec& m)
{
  double val=0;
  cpx* data = m.data();
  for(int i=0; i<m.m(); i++) {	 //val += data[i].re*data[i].re + data[i].im*data[i].im;
	 val += norm(data[i]);
  }
  return val;
}
inline double energy(DblOffVec& m)
{
  double val=0;
  double* data = m.data();
  for(int i=0; i<m.m(); i++) {
	 val += data[i]*data[i];
  }
  return val;
}
inline double energy(DblNumVec& m)
{
  double val=0;
  double* data = m.data();
  for(int i=0; i<m.m(); i++) {
	 val += data[i]*data[i];
  }
  return val;
}

FDCT_USFFT_NS_END_NAMESPACE

#endif
