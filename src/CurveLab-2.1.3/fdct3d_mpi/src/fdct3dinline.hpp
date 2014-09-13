/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#ifndef _FDCT3DINLINE_HPP_
#define _FDCT3DINLINE_HPP_

#include "fdct3d.hpp"

//-------------------------------------------------
inline int fdct3d_position_aux(int N1,int N2,int N3,int b, int x,int y,int z,
									 int& bi,int& bj,int& bk, int& oi,int& oj,int& ok)
{
  x += N1/2;  y += N2/2;  z += N3/2;
  bi = x/b;  bj = y/b;  bk = z/b;
  oi = x%b;  oj = y%b;  ok = z%b;
  return 0;
}

//-------------------------------------------------
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

//-------------------------------------------------
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

//-------------------------------------------------
inline int fdct3d_rangecompute(double L1, double L2, double L3,
									 int& S1, int& S2, int& S3,
									 int& F1, int& F2, int& F3,
									 double& R1, double& R2, double& R3)
{
  S1 = 2*int(floor(L1/2))+1;  S2 = 2*int(floor(L2/2))+1;  S3 = 2*int(floor(L3/2))+1; //number of samples
  F1 = int(floor(L1/2));  F2 = int(floor(L2/2));  F3 = int(floor(L3/2)); //offset on either side
  R1 = L1/2;  R2 = L2/2;  R3 = L3/2;
  return 0;
}

//-------------------------------------------------
//windowing
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

#endif

