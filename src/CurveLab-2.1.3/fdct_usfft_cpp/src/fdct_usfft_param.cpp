/*
  Copyright (C) 2004 Caltech
  Written by Lexing Ying
*/

#include "fdct_usfft.hpp"
#include "fdct_usfft_inline.hpp"

FDCT_USFFT_NS_BEGIN_NAMESPACE

int fdct_usfft_param_sepangle(double XL1, double XL2, int nbangle,
										vector<dbl2>& sx, vector<dbl2>& sy,
										vector<double>& fx, vector<double>& fy,
										vector<int>& nx, vector<int>& ny);
int fdct_usfft_param_wavelet(int N1, int N2,
									  vector<dbl2>& sx, vector<dbl2>& sy,
									  vector<double>& fx, vector<double>& fy,
									  vector<int>& nx, vector<int>& ny);

//----------------------------------------------------------------------
int fdct_usfft_param(int N1, int N2, int nbscales, int nbangles_coarse, int ac,
							vector< vector<dbl2> >& sx, vector< vector<dbl2> >& sy,
							vector< vector<double> >& fx, vector< vector<double> >& fy,
							vector< vector<int> >& nx, vector< vector<int> >& ny)
{
  sx.resize(nbscales);  sy.resize(nbscales);
  fx.resize(nbscales);  fy.resize(nbscales);
  nx.resize(nbscales);  ny.resize(nbscales);
  
  if(ac==1) {
	 //nbangles
	 vector<int> nbangles(nbscales);
	 nbangles[0] = 1;
	 for(int sc=1; sc<nbscales; sc++)		nbangles[sc] = nbangles_coarse * pow2( int(ceil(double(sc-1)/2)) );
	 
	 //high freq levels
	 double XL1 = 4.0*N1/3.0;  double XL2 = 4.0*N2/3.0; //range
	 for(int sc=nbscales-1; sc>0; sc--) {
		fdct_usfft_param_sepangle(XL1, XL2, nbangles[sc], sx[sc], sy[sc], fx[sc], fy[sc], nx[sc], ny[sc]);
		XL1 /= 2;		XL2 /= 2;
	 }
	 //coarsest level
	 int XS1, XS2;  int XF1, XF2;  double XR1, XR2;	 fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
	 fdct_usfft_param_wavelet(XS1, XS2, sx[0], sy[0], fx[0], fy[0], nx[0], ny[0]);
  } else {
	 //nbangles
	 vector<int> nbangles(nbscales);
	 nbangles[0] = 1;
	 for(int sc=1; sc<nbscales-1; sc++)		nbangles[sc] = nbangles_coarse * pow2( int(ceil(double(sc-1)/2)) );
	 nbangles[nbscales-1] = 1;
	 //top level
	 fdct_usfft_param_wavelet(N1, N2, sx[nbscales-1], sy[nbscales-1], fx[nbscales-1], fy[nbscales-1], nx[nbscales-1], ny[nbscales-1]);
	 //next levels
	 double XL1 = 2.0*N1/3.0;  double XL2 = 2.0*N2/3.0; //range
	 for(int sc=nbscales-2; sc>0; sc--) {
		fdct_usfft_param_sepangle(XL1, XL2, nbangles[sc], sx[sc], sy[sc], fx[sc], fy[sc], nx[sc], ny[sc]);
		XL1 /= 2;		XL2 /= 2;
	 }
	 //coarsest level
	 int XS1, XS2;  int XF1, XF2;  double XR1, XR2;	 fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
	 fdct_usfft_param_wavelet(XS1, XS2, sx[0], sy[0], fx[0], fy[0], nx[0], ny[0]);
  }
  
  return 0;
}

int fdct_usfft_param_sepangle(double XL1, double XL2, int nbangle,
										vector<dbl2>& sx, vector<dbl2>& sy,
										vector<double>& fx, vector<double>& fy,
										vector<int>& nx, vector<int>& ny)
{
  fx.resize(nbangle);  fy.resize(nbangle);
  nx.resize(nbangle);  ny.resize(nbangle);
  sx.resize(nbangle);  sy.resize(nbangle);
  
  int nd = nbangle / 4;
  int wcnt = 0;
  
  int XS1, XS2;  int XF1, XF2;  double XR1, XR2;	 fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
  double XW1 = XL1/nd;	 double XW2 = XL2/nd;
  //2.
  for(int w=nd-1; w>=0; w--) {
	 double xs = -XR1;		double xe = -XR1/4;		int xn = int(ceil(xe-xs));		int yn = 2*int(ceil(XW2))+1;
	 double ym = XR2 - (w+0.5)*XW2;
	 double s2 = -ym/XR1;
	 if(xn%2==0) xn++;	 if(yn%2==0) yn++;
	 fx[wcnt] = -XR1/2;	 fy[wcnt] = ym/2;
	 nx[wcnt] = xn;	 ny[wcnt] = yn;
	 sx[wcnt] = dbl2(1.0/xn,0);	 sy[wcnt] = dbl2(-s2*1.0/yn, 1.0/yn);
	 wcnt++;
  }
  //1.
  for(int w=nd-1; w>=0; w--) {
	 double ys = XR2/4;	 double ye = XR2;	 int yn = int(ceil(ye-ys));	 int xn = 2*int(ceil(XW1))+1;
	 double xm = XR1 - (w+0.5)*XW1;
	 double s2 = -xm/XR2;
	 if(xn%2==0) xn++;	 if(yn%2==0) yn++;
	 fx[wcnt] = xm/2;	 fy[wcnt] = XR2/2;
	 nx[wcnt] = xn;	 ny[wcnt] = yn;
	 sx[wcnt] = dbl2(1.0/xn, s2*1.0/xn);	 sy[wcnt] = dbl2(0, 1.0/yn);
	 wcnt++;
  }
  //0.
  for(int w=nd-1; w>=0; w--) {
  	 double xs = XR1/4;	 double xe = XR1;	 int xn = int(ceil(xe-xs));	 int yn = 2*int(ceil(XW2))+1;
	 double ym = -XR2 + (w+0.5)*XW2;
	 double s2 = ym/XR1;
	 if(xn%2==0) xn++;	 if(yn%2==0) yn++;
	 fx[wcnt] = XR1/2;	 fy[wcnt] = ym/2;
	 nx[wcnt] = xn;	 ny[wcnt] = yn;
	 sx[wcnt] = dbl2(1.0/xn,0);	 sy[wcnt] = dbl2(-s2*1.0/yn, 1.0/yn);
	 wcnt++;
  }
  //3.
  for(int w=nd-1; w>=0; w--) {
	 double ys = -XR2;	 double ye = -XR2/4;	 int yn = int(ceil(ye-ys));	 int xn = 2*int(ceil(XW1))+1;
	 double xm = -XR1 + (w+0.5)*XW1;
	 double s2 = xm/XR2;
	 if(xn%2==0) xn++;	 if(yn%2==0) yn++;
	 fx[wcnt] = xm/2;	 fy[wcnt] = -XR2/2;
	 nx[wcnt] = xn;	 ny[wcnt] = yn;
	 sx[wcnt] = dbl2(1.0/xn, s2*1.0/xn);	 sy[wcnt] = dbl2(0, 1.0/yn);
	 wcnt++;
  }
  assert(wcnt==nbangle);
  return 0;
}

int fdct_usfft_param_wavelet(int N1, int N2,
									  vector<dbl2>& sx, vector<dbl2>& sy,
									  vector<double>& fx, vector<double>& fy,
									  vector<int>& nx, vector<int>& ny)
{
  fx.resize(1);  fx[0] = 0;
  fy.resize(1);  fy[0] = 0;
  nx.resize(1);  nx[0] = N1;
  ny.resize(1);  ny[0] = N2;
  double dx = 1.0/N1;  double dy = 1.0/N2;
  sx.resize(1);  sx[0] = dbl2(dx,0);
  sy.resize(1);  sy[0] = dbl2(0,dy);
  return 0;
}

FDCT_USFFT_NS_END_NAMESPACE
