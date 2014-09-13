/*
   Copyright (C) 2004 Caltech
   Written by Lexing Ying
*/

#include "fdct_wrapping.hpp"
#include "fdct_wrapping_inline.hpp"

FDCT_WRAPPING_NS_BEGIN_NAMESPACE

int fdct_wrapping_param_sepangle(double XL1, double XL2, int nbangle,
											vector<double>& sx, vector<double>& sy,
											vector<double>& fx, vector<double>& fy,
											vector<int>& nx, vector<int>& ny);
int fdct_wrapping_param_wavelet(int N1, int N2,
										  vector<double>& sx, vector<double>& sy,
										  vector<double>& fx, vector<double>& fy,
										  vector<int>& nx, vector<int>& ny);

//--------------------------------------------------
int fdct_wrapping_param(int N1, int N2, int nbscales, int nbangles_coarse, int allcurvelets,
								vector< vector<double> >& sx, vector< vector<double> >& sy,
								vector< vector<double> >& fx, vector< vector<double> >& fy,
								vector< vector<int> >& nx, vector< vector<int> >& ny)
{
  //sx, sy, step in spatial domain
  //fx, fy, position in frequency domain
  //nx, ny, size of the grid
  sx.resize(nbscales);  sy.resize(nbscales);
  fx.resize(nbscales);  fy.resize(nbscales);
  nx.resize(nbscales);  ny.resize(nbscales);
    
  vector<int> nbangles(nbscales);
  if(allcurvelets==1) {
	 //nbangles
	 nbangles[0] = 1;
	 for(int sc=1; sc<nbscales; sc++)		nbangles[sc] = nbangles_coarse * pow2( int(ceil(double(sc-1)/2)) );
	 //high freq levels
	 double XL1 = 4.0*N1/3.0;  double XL2 = 4.0*N2/3.0; //range
	 for(int sc=nbscales-1; sc>0; sc--) {
		fdct_wrapping_param_sepangle(XL1, XL2, nbangles[sc], sx[sc], sy[sc], fx[sc], fy[sc], nx[sc], ny[sc]);
		XL1 /= 2;		XL2 /= 2;
	 }
	 //coarsest level
	 int XS1, XS2;  int XF1, XF2;  double XR1, XR2;	 fdct_wrapping_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
	 fdct_wrapping_param_wavelet(XS1, XS2, sx[0], sy[0], fx[0], fy[0], nx[0], ny[0]);
  } else {
	 //nbangles
	 nbangles[0] = 1;
	 for(int sc=1; sc<nbscales-1; sc++)		nbangles[sc] = nbangles_coarse * pow2( int(ceil(double(sc-1)/2)) );
	 nbangles[nbscales-1] = 1;
	 //top level
	 fdct_wrapping_param_wavelet(N1, N2, sx[nbscales-1], sy[nbscales-1], fx[nbscales-1], fy[nbscales-1], nx[nbscales-1], ny[nbscales-1]);
	 //next levels
	 double XL1 = 2.0*N1/3.0;  double XL2 = 2.0*N2/3.0; //range
	 for(int sc=nbscales-2; sc>0; sc--) {
		fdct_wrapping_param_sepangle(XL1, XL2, nbangles[sc], sx[sc], sy[sc], fx[sc], fy[sc], nx[sc], ny[sc]);
		XL1 /= 2;		XL2 /= 2;
	 }
	 //coarsest level
	 int XS1, XS2;  int XF1, XF2;  double XR1, XR2;	 fdct_wrapping_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
	 fdct_wrapping_param_wavelet(XS1, XS2, sx[0], sy[0], fx[0], fy[0], nx[0], ny[0]);
  }
  
  return 0;
}

//--------------------------------------------------
int fdct_wrapping_param_sepangle(double XL1, double XL2, int nbangle,
											vector<double>& sx, vector<double>& sy,
											vector<double>& fx, vector<double>& fy,
											vector<int>& nx, vector<int>& ny)
{
  fx.resize(nbangle);  fy.resize(nbangle);
  nx.resize(nbangle);  ny.resize(nbangle);
  sx.resize(nbangle);  sy.resize(nbangle);
 
  int nbquadrants = 4;
  int nd = nbangle / 4;  //int nbangles_perquad = nbangles[sc] / 4;
  int wcnt = 0;
  //backup
  double XL1b = XL1;  double XL2b = XL2;
  int qvec[] = {2,1,0,3};
  for(int qi=0; qi<nbquadrants; qi++) {
	 int q = qvec[qi];
	 fdct_wrapping_rotate_forward(q, XL1b, XL2b, XL1, XL2);	 XL1 = abs(XL1);	 XL2 = abs(XL2);
	 double XW1 = XL1/nd;	 double XW2 = XL2/nd;
	 int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_wrapping_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
	 for(int w=nd-1; w>=0; w--) {
		//get size
		double xs = XR1/4 - (XW1/2)/4;
		double xe = XR1;
		double ys = -XR2 + (w-0.5)*XW2;
		double ye = -XR2 + (w+1.5)*XW2; //x range
		int xn = int(ceil(xe-xs));			 int yn = int(ceil(ye-ys));
		//MAKE THEM ODD
		if(xn%2==0) xn++;		if(yn%2==0) yn++;
		
		double tx; double ty;
		//fx,fy
		tx = XR1/2;		ty = (-XR2 + (w+0.5)*XW2)/2; //center, NOTE: need /2
		fdct_wrapping_rotate_backward(q, tx, ty, fx[wcnt], fy[wcnt]);
		//sx, sy;
		tx = 1.0/xn;		ty = 1.0/yn;
		fdct_wrapping_rotate_backward(q, xn, yn, nx[wcnt], ny[wcnt]);		nx[wcnt] = abs(nx[wcnt]);		ny[wcnt] = abs(ny[wcnt]);
		fdct_wrapping_rotate_backward(q, tx, ty, sx[wcnt], sy[wcnt]);		sx[wcnt] = abs(sx[wcnt]);		sy[wcnt] = abs(sy[wcnt]);
		wcnt++;
	 }
  }
  
  return 0;
}

//--------------------------------------------------
int fdct_wrapping_param_wavelet(int S1, int S2,
										  vector<double>& sx, vector<double>& sy,
										  vector<double>& fx, vector<double>& fy,
										  vector<int>& nx, vector<int>& ny)
{
  fx.resize(1);  fx[0] = 0;
  fy.resize(1);  fy[0] = 0;
  nx.resize(1);  nx[0] = S1;
  ny.resize(1);  ny[0] = S2;
  double dx = 1.0/S1;  double dy = 1.0/S2;
  sx.resize(1);  sx[0] = dx;
  sy.resize(1);  sy[0] = dy;
  return 0;
}

FDCT_WRAPPING_NS_END_NAMESPACE
