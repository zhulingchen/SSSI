/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#ifndef _FDCT3DINLINE_HPP_
#define _FDCT3DINLINE_HPP_

#include "fdct3d.hpp"

inline int fdct3d_fftshift(int N1, int N2, int N3, CpxNumTns& T, CpxOffTns& O)
{
  int F1 = N1/2;  int F2 = N2/2;   int F3 = N3/2; 
  IntOffVec p1(N1);  for(int i=-F1; i<-F1+N1; i++)	 p1(i) = (i+N1)%N1;
  IntOffVec p2(N2);  for(int i=-F2; i<-F2+N2; i++)	 p2(i) = (i+N2)%N2;
  IntOffVec p3(N3);  for(int i=-F3; i<-F3+N3; i++)	 p3(i) = (i+N3)%N3;
  for(int i=-F1; i<-F1+N1; i++)
	 for(int j=-F2; j<-F2+N2; j++)
		for(int k=-F3; k<-F3+N3; k++)
		  O(i,j,k) = T(p1(i),p2(j),p3(k));
  return 0;
}

inline int fdct3d_ifftshift(int N1, int N2, int N3, CpxOffTns& O, CpxNumTns& T)
{
  int F1 = N1/2;  int F2 = N2/2;   int F3 = N3/2; 
  IntOffVec p1(N1);  for(int i=-F1; i<-F1+N1; i++)	 p1(i) = (i+N1)%N1;
  IntOffVec p2(N2);  for(int i=-F2; i<-F2+N2; i++)	 p2(i) = (i+N2)%N2;
  IntOffVec p3(N3);  for(int i=-F3; i<-F3+N3; i++)	 p3(i) = (i+N3)%N3;
  for(int i=-F1; i<-F1+N1; i++)
	 for(int j=-F2; j<-F2+N2; j++)
		for(int k=-F3; k<-F3+N3; k++)
		  T(p1(i),p2(j),p3(k)) = O(i,j,k);
  return 0;
}

inline int fdct3d_window(double x, double& l, double& r)
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

inline int fdct3d_rangecompute(double XL1, double XL2, double XL3,
									 int& XS1, int& XS2, int& XS3,
									 int& XF1, int& XF2, int& XF3,
									 double& XR1, double& XR2, double& XR3)
{
  XS1 = 2*int(floor(XL1/2))+1;
  XS2 = 2*int(floor(XL2/2))+1;
  XS3 = 2*int(floor(XL3/2))+1; //number of samples
  
  XF1 = int(floor(XL1/2));
  XF2 = int(floor(XL2/2));
  XF3 = int(floor(XL3/2)); //offset on either side
  
  XR1 = XL1/2;
  XR2 = XL2/2;
  XR3 = XL3/2;
  return 0;
}

inline int fdct3d_lowpass(double L1, DblOffVec& lowpass)
{
  int S1 = 2*int(floor(L1/2))+1;
  int F1 = int(floor(L1/2));
  double R1 = L1/2;
  
  assert(lowpass.m()>=S1);
  setvalue(lowpass, double(0.0));
  for(int i=-F1; i<=0; i++) {
	 double x = (i+R1)/(R1/2);
	 double l,r; fdct3d_window(x, l, r);
	 lowpass(i) = l;
	 lowpass(-i) = l;
  }
  
  return 0;
}

inline int fdct3d_lowpasscompute(double XL1, double XL2, double XL3, DblOffTns& lowpass)
{
  int XS1 = 2*int(floor(XL1/2))+1;
  int XS2 = 2*int(floor(XL2/2))+1;
  int XS3 = 2*int(floor(XL3/2))+1;
  //number of samples
  int XF1 = int(floor(XL1/2));
  int XF2 = int(floor(XL2/2));
  int XF3 = int(floor(XL3/2));
  //offset on either side
  double XR1 = XL1/2;
  double XR2 = XL2/2;
  double XR3 = XL3/2;
  
  DblOffVec lowpass1(XS1); setvalue(lowpass1,1.0);
  for(int i=-XF1; i<-XR1/2; i++) {
	 double x = (i+XR1)/(XR1/2);
	 double l,r; fdct3d_window(x, l, r);
	 lowpass1(i) = l;
	 lowpass1(-i) = l;
  }  //cerr<<lowpass1;
  DblOffVec lowpass2(XS2); setvalue(lowpass2,1.0);
  for(int i=-XF2; i<-XR2/2; i++) {
	 double x = (i+XR2)/(XR2/2);
	 double l,r; fdct3d_window(x, l, r);
	 lowpass2(i) = l;
	 lowpass2(-i) = l;
  }  //cerr<<lowpass2;
  DblOffVec lowpass3(XS3); setvalue(lowpass3,1.0);
  for(int i=-XF3; i<-XR3/2; i++) {
	 double x = (i+XR3)/(XR3/2);
	 double l,r; fdct3d_window(x, l, r);
	 lowpass3(i) = l;
	 lowpass3(-i) = l;
  }  //cerr<<lowpass3;
  
  for(int i=-XF1; i<-XF1+XS1; i++)
	 for(int j=-XF2; j<-XF2+XS2; j++)
		for(int k=-XF3; k<-XF3+XS3; k++)
		  lowpass(i,j,k) = lowpass1(i) * lowpass2(j) * lowpass3(k);
  
  return 0;
}

template<class F>
inline int fdct3d_rotate_forward(int f, F XL1, F XL2, F XL3, F& TL1, F& TL2, F& TL3)
{
  if(f==0) {
	 TL1 = XL1;	 TL2 = XL2;  TL3 = XL3;
  } else if(f==1) {
	 TL1 = -XL1;  TL2 = -XL2;  TL3 = -XL3;
  } else if(f==2) {
	 TL1 = XL3; TL2 = XL1; TL3 = XL2;
  } else if(f==3) {
	 TL1 = -XL3; TL2 = -XL1; TL3 = -XL2;
  } else if(f==4) {
	 TL1 = XL2; TL2 = XL3; TL3 = XL1;
  } else if(f==5) {
	 TL1 = -XL2; TL2 = -XL3; TL3 = -XL1;
  }
  return 0;
}

template<class F>
inline int fdct3d_rotate_backward(int f, F XL1, F XL2, F XL3, F& TL1, F& TL2, F& TL3)
{
  if(f==0) {
	 TL1 = XL1;	 TL2 = XL2;  TL3 = XL3;
  } else if(f==1) {
	 TL1 = -XL1;  TL2 = -XL2;  TL3 = -XL3;
  } else if(f==2) {
	 TL1 = XL2; TL2 = XL3; TL3 = XL1;
  } else if(f==3) {
	 TL1 = -XL2; TL2 = -XL3; TL3 = -XL1;
  } else if(f==4) {
	 TL1 = XL3; TL2 = XL1; TL3 = XL2;
  } else if(f==5) {
	 TL1 = -XL3; TL2 = -XL1; TL3 = -XL2;
  }
  return 0;
}

inline int fdct3d_rotate_forward(int f, CpxOffTns& X, CpxOffTns& T)
{
  //dihedral group
  int XS1 = X.m();  int XS2 = X.n();  int XS3 = X.p();
  int XF1 =-X.s();  int XF2 =-X.t();  int XF3 =-X.u();
  int TS1, TS2, TS3;  //int TF1, TF2, TF3;
  if(f==0) {
	 TS1 = XS1;	 TS2 = XS2;	 TS3 = XS3;
	 T.resize(TS1, TS2, TS3);
	 for(int i=-XF1; i<-XF1+XS1; i++) {
		for(int j=-XF2; j<-XF2+XS2; j++) {
		  for(int k=-XF3; k<-XF3+XS3; k++) {
			 T(i,j,k) = X(i,j,k);
		  }
		}
	 }
  } else if(f==1) {
	 TS1 = XS1;	 TS2 = XS2;	 TS3 = XS3;
	 T.resize(TS1, TS2, TS3);
	 for(int i=-XF1; i<-XF1+XS1; i++) {		//int ni = (-i<-XF1+XS1) ? -i : -i-XS1;
		for(int j=-XF2; j<-XF2+XS2; j++) {		  //int nj = (-j<-XF2+XS2) ? -j : -j-XS2;
		  for(int k=-XF3; k<-XF3+XS3; k++) {			 //int nk = (-k<-XF3+XS3) ? -k : -k-XS3;
			 T(-i,-j,-k) = X(i,j,k);
		  }
		}
	 }
  } else if(f==2) {
	 TS1 = XS3; TS2 = XS1; TS3 = XS2;
	 T.resize(TS1, TS2, TS3);
	 for(int i=-XF1; i<-XF1+XS1; i++) {
		for(int j=-XF2; j<-XF2+XS2; j++) {
		  for(int k=-XF3; k<-XF3+XS3; k++) {
			 T(k,i,j) = X(i,j,k);
		  }
		}
	 }
  } else if(f==3) {
	 TS1 = XS3; TS2 = XS1; TS3 = XS2;
	 T.resize(TS1, TS2, TS3);
	 for(int i=-XF1; i<-XF1+XS1; i++) {		//int ni = (-i<-XF1+XS1) ? -i : -i-XS1;
		for(int j=-XF2; j<-XF2+XS2; j++) {		  //int nj = (-j<-XF2+XS2) ? -j : -j-XS2;
		  for(int k=-XF3; k<-XF3+XS3; k++) {			 //int nk = (-k<-XF3+XS3) ? -k : -k-XS3;
			 T(-k,-i,-j) = X(i,j,k);
		  }
		}
	 }
  } else if(f==4) {
	 TS1 = XS2; TS2 = XS3; TS3 = XS1;	 //XF1 = XF2; XF2 = XF3; XF3 = XF1;	 //TL1 = XL2; TL2 = XL3; TL3 = XL1;
	 T.resize(TS1, TS2, TS3);
	 for(int i=-XF1; i<-XF1+XS1; i++) {
		for(int j=-XF2; j<-XF2+XS2; j++) {
		  for(int k=-XF3; k<-XF3+XS3; k++) {
			 T(j,k,i) = X(i,j,k);
		  }
		}
	 }
  } else if(f==5) {
	 TS1 = XS2; TS2 = XS3; TS3 = XS1;	 //XF1 = XF3; XF2 = XF2; XF3 = XF1;	 //TL1 = XL3; TL2 = XL2; TL3 = XL1;
	 T.resize(TS1, TS2, TS3);
	 for(int i=-XF1; i<-XF1+XS1; i++) {		//int ni = (-i<-XF1+XS1) ? -i : -i-XS1;
		for(int j=-XF2; j<-XF2+XS2; j++) {		  //int nj = (-j<-XF2+XS2) ? -j : -j-XS2;
		  for(int k=-XF3; k<-XF3+XS3; k++) {			 //int nk = (-k<-XF3+XS3) ? -k : -k-XS3;
			 T(-j,-k,-i) = X(i,j,k);
		  }
		}
	 }
  }
  return 0;
}

inline int fdct3d_rotate_backward(int f, CpxOffTns& X, CpxOffTns& T)
{
  //dihedral group
  int XS1 = X.m();  int XS2 = X.n();  int XS3 = X.p();
  int XF1 =-X.s();  int XF2 =-X.t();  int XF3 =-X.u();
  int TS1, TS2, TS3;  //int XF1, XF2, XF3;
  if(f==0) {
	 TS1 = XS1;	 TS2 = XS2;	 TS3 = XS3;	 //XF1 = XF1;	 XF2 = XF2;	 XF3 = XF3;	 //TL1 = XL1;	 TL2 = XL2;  TL3 = XL3;
	 T.resize(TS1, TS2, TS3);
	 for(int i=-XF1; i<-XF1+XS1; i++) {
		for(int j=-XF2; j<-XF2+XS2; j++) {
		  for(int k=-XF3; k<-XF3+XS3; k++) {
			 T(i,j,k) = X(i,j,k);
		  }
		}
	 }
  } else if(f==1) {
	 TS1 = XS1;	 TS2 = XS2;	 TS3 = XS3;
	 T.resize(TS1, TS2, TS3);
	 for(int i=-XF1; i<-XF1+XS1; i++) {		//int ni = (-i<-XF1+XS1) ? -i : -i-TS1;
		for(int j=-XF2; j<-XF2+XS2; j++) {		  //int nj = (-j<-XF2+XS2) ? -j : -j-TS2;
		  for(int k=-XF3; k<-XF3+XS3; k++) {			 //int nk = (-k<-XF3+XS3) ? -k : -k-TS3;
			 T(-i,-j,-k) = X(i,j,k);
		  }
		}
	 }
  } else if(f==2) {
	 TS1 = XS2; TS2 = XS3; TS3 = XS1;	 //XF1 = XF2; XF2 = XF3; XF3 = XF1;	 //TL1 = XL2; TL2 = XL3; TL3 = XL1;
	 T.resize(TS1, TS2, TS3);
	 for(int i=-XF1; i<-XF1+XS1; i++) {
		for(int j=-XF2; j<-XF2+XS2; j++) {
		  for(int k=-XF3; k<-XF3+XS3; k++) {
			 T(j,k,i) = X(i,j,k);
		  }
		}
	 }
  } else if(f==3) {
	 TS1 = XS2; TS2 = XS3; TS3 = XS1;	 //XF1 = XF2; XF2 = XF1; XF3 = XF3;	 //TL1 = XL2; TL2 = XL1; TL3 = XL3;
	 T.resize(TS1, TS2, TS3);
	 for(int i=-XF1; i<-XF1+XS1; i++) {		//int ni = (-i<-XF1+XS1) ? -i : -i-TS1;
		for(int j=-XF2; j<-XF2+XS2; j++) {		  //int nj = (-j<-XF2+XS2) ? -j : -j-TS2;
		  for(int k=-XF3; k<-XF3+XS3; k++) {			 //int nk = (-k<-XF3+XS3) ? -k : -k-TS3;
			 T(-j,-k,-i) = X(i,j,k);
		  }
		}
	 }
  } else if(f==4) {
	 TS1 = XS3; TS2 = XS1; TS3 = XS2;	 //XF1 = XF3; XF2 = XF1; XF3 = XF2;	 //TL1 = XL3; TL2 = XL1; TL3 = XL2;
	 T.resize(TS1, TS2, TS3);
	 for(int i=-XF1; i<-XF1+XS1; i++) {
		for(int j=-XF2; j<-XF2+XS2; j++) {
		  for(int k=-XF3; k<-XF3+XS3; k++) {
			 T(k,i,j) = X(i,j,k);
		  }
		}
	 }
  } else if(f==5) {
	 TS1 = XS3; TS2 = XS1; TS3 = XS2;	 //XF1 = XF3; XF2 = XF2; XF3 = XF1;	 //TL1 = XL3; TL2 = XL2; TL3 = XL1;
	 T.resize(TS1, TS2, TS3);
	 for(int i=-XF1; i<-XF1+XS1; i++) {		//int ni = (-i<-XF1+XS1) ? -i : -i-TS1;
		for(int j=-XF2; j<-XF2+XS2; j++) {		  //int nj = (-j<-XF2+XS2) ? -j : -j-TS2;
		  for(int k=-XF3; k<-XF3+XS3; k++) {			 //int nk = (-k<-XF3+XS3) ? -k : -k-TS3;
			 T(-k,-i,-j) = X(i,j,k);
		  }
		}
	 }
  }
  return 0;
}

inline int fdct3d_globalpou_shape(double lt, double rt, double u, double v, double& s)
{
  double t,a,b;
  t = (u-lt)/(rt-lt);
  if(t<=0) a = 1.0; else if(t>=1) a = 0.0; else a = exp(2.0*exp(-1.0/(t))/(t-1.0));
  t = (v-lt)/(rt-lt);
  if(t<=0) b = 1.0; else if(t>=1) b = 0.0; else b = exp(2.0*exp(-1.0/(t))/(t-1.0));
  s = a*b;
  return 0;
}

inline int fdct3d_globalpou(double tht, double phi, double buf, double& pou)
{
  double lt = M_PI/4 - buf;  double rt = M_PI/4 + buf;
  if(abs(tht)<=lt && abs(phi)<=lt)
	 pou = 1;
  else if(abs(tht)>=rt || abs(phi)>=rt)
	 pou = 0;
  else {
	 double tt = tan(abs(tht));
	 double tp = tan(abs(phi));
	 double pa[2];	pa[0] = abs(tht);	pa[1] = abs(phi);
	 double pb[2]; pb[0] = atan2(tp,tt); pb[1] = atan2(1,tt);
	 double pc[2]; pc[0] = atan2(tt,tp); pc[1] = atan2(1,tp);
	 double sa; fdct3d_globalpou_shape(lt, rt, pa[0], pa[1], sa);
	 double sb; fdct3d_globalpou_shape(lt, rt, pb[0], pb[1], sb);
	 double sc; fdct3d_globalpou_shape(lt, rt, pc[0], pc[1], sc);
	 pou = sqrt(sa/(sa+sb+sc));
  }
  return 0;
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
inline double energy(CpxOffMat& m)
{
  double val=0;
  cpx* data = m.data();
  for(int i=0; i<m.m()*m.n(); i++) {	 //val += data[i].re*data[i].re + data[i].im*data[i].im;
	 val += norm(data[i]);
  }
  return val;
}
inline double energy(CpxOffTns& m)
{
  double val=0;
  cpx* data = m.data();
  for(int i=0; i<m.m()*m.n()*m.p(); i++) {	 //val += data[i].re*data[i].re + data[i].im*data[i].im;
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
inline double energy(CpxNumMat& m)
{
  double val=0;
  cpx* data = m.data();
  for(int i=0; i<m.m()*m.n(); i++) {	 //val += data[i].re*data[i].re + data[i].im*data[i].im;
	 val += norm(data[i]);
  }
  return val;
}
inline double energy(CpxNumTns& m)
{
  double val=0;
  cpx* data = m.data();
  for(int i=0; i<m.m()*m.n()*m.p(); i++) {	 //val += data[i].re*data[i].re + data[i].im*data[i].im;
	 val += norm(data[i]);
  }
  return val;
}


#endif
