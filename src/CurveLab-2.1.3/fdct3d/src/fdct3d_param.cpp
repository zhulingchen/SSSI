/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/
#include "fdct3d.hpp"
#include "fdct3dinline.hpp"

//----------------------------------------
int fdct3d_param_center(double L1, double L2, double L3, int s,
							vector< vector<double> >& fxs, vector< vector<double> >& fys, vector< vector<double> >& fzs,
							vector< vector<int   > >& nxs, vector< vector<int   > >& nys, vector< vector<int   > >& nzs)
{
  vector<double>& fx = fxs[s];  vector<double>& fy = fys[s];  vector<double>& fz = fzs[s];
  vector<int>& nx = nxs[s];  vector<int>& ny = nys[s];  vector<int>& nz = nzs[s];
  
  int S1, S2, S3;	 int F1, F2, F3;	 double R1, R2, R3;	 fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);
  fx.resize(1);  fy.resize(1);  fz.resize(1);
  nx.resize(1);  ny.resize(1);  nz.resize(1);
  
  fx[0] = 0;  fy[0] = 0;  fz[0] = 0;
  nx[0] = S1;  ny[0] = S2;  nz[0] = S3;
  
  return 0;
}

int fdct3d_param_wavelet(double L1, double L2, double L3, int s, int N1,int N2,int N3,
							 vector< vector<double> >& fxs, vector< vector<double> >& fys, vector< vector<double> >& fzs,
							 vector< vector<int   > >& nxs, vector< vector<int   > >& nys, vector< vector<int   > >& nzs)
{
  vector<double>& fx = fxs[s];  vector<double>& fy = fys[s];  vector<double>& fz = fzs[s];
  vector<int>& nx = nxs[s];  vector<int>& ny = nys[s];  vector<int>& nz = nzs[s];
  
  fx.resize(1);  fy.resize(1);  fz.resize(1);
  nx.resize(1);  ny.resize(1);  nz.resize(1);
  
  fx[0] = 0;  fy[0] = 0;  fz[0] = 0;
  nx[0] = N1;  ny[0] = N2;  nz[0] = N3;

  return 0;
}

//----------------------------------------
int fdct3d_param_angles(double L1, double L2, double L3, int s, int nd,
							vector< vector<double> >& fxs, vector< vector<double> >& fys, vector< vector<double> >& fzs,
							vector< vector<int   > >& nxs, vector< vector<int   > >& nys, vector< vector<int   > >& nzs)
{
  vector<double>& fx = fxs[s];  vector<double>& fy = fys[s];  vector<double>& fz = fzs[s];
  vector<int>& nx = nxs[s];  vector<int>& ny = nys[s];  vector<int>& nz = nzs[s];
  
  int nbw = 6 * nd * nd;
  fx.resize(nbw);  fy.resize(nbw);  fz.resize(nbw);
  nx.resize(nbw);  ny.resize(nbw);  nz.resize(nbw);
  
  int nf = 6;
  int wcnt = 0;
  int S1, S2, S3;	 int F1, F2, F3;	 double R1, R2, R3;	 fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);
  double W1 = L1/nd;  double W2 = L2/nd;  double W3 = L3/nd;

  //face 0: x,y,z
  for(int h=0; h<nd; h++) { //(y first z second)
	 for(int g=0; g<nd; g++) {
		double xs = R1/4-(W1/2)/4;		double xe = R1;
		double ys = -R2 + (2*g-1)*W2/2;		double ye = -R2 + (2*g+3)*W2/2;
		double zs = -R3 + (2*h-1)*W3/2;		double ze = -R3 + (2*h+3)*W3/2;
		int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		double tx, ty, tz;
		tx = R1/2;		  ty = (-R2 + (g+0.5)*W2)/2;		  tz = (-R3 + (h+0.5)*W3)/2;
		fx[wcnt] = tx;		  fy[wcnt] = ty;		  fz[wcnt] = tz;
		nx[wcnt] = xn;		  ny[wcnt] = yn;		  nz[wcnt] = zn;
		wcnt++;
	 }
  }
  //face 1: y z x
  for(int f=0; f<nd; f++) {
	 for(int h=0; h<nd; h++) {
		double ys = R2/4-(W2/2)/4;		  double ye = R2;
		double zs = -R3 + (2*h-1)*W3/2;		  double ze = -R3 + (2*h+3)*W3/2;
		double xs = -R1 + (2*f-1)*W1/2;		  double xe = -R1 + (2*f+3)*W1/2;
		int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		double tx, ty, tz;
		ty = R2/2;		  tz = (-R3 + (h+0.5)*W3)/2;		  tx = (-R1 + (f+0.5)*W1)/2;
		fx[wcnt] = tx;		  fy[wcnt] = ty;		  fz[wcnt] = tz;
		nx[wcnt] = xn;		  ny[wcnt] = yn;		  nz[wcnt] = zn;
		wcnt++;
	 }
  }
  //face 2: z,x,y
  for(int g=0; g<nd; g++) {
	 for(int f=0; f<nd; f++) {
		double zs = R3/4-(W3/2)/4;		double ze = R3;
		double xs = -R1 + (2*f-1)*W1/2;		double xe = -R1 + (2*f+3)*W1/2;
		double ys = -R2 + (2*g-1)*W2/2;		double ye = -R2 + (2*g+3)*W2/2;
		int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		double tx, ty, tz;
		tz = R3/2;		  tx = (-R1 + (f+0.5)*W1)/2;		  ty = (-R2 + (g+0.5)*W2)/2;
		fx[wcnt] = tx;		  fy[wcnt] = ty;		  fz[wcnt] = tz;
		nx[wcnt] = xn;		  ny[wcnt] = yn;		  nz[wcnt] = zn;
		wcnt++;
	 }
  }
  //face 3: -x,-y,-z
  for(int h=nd-1; h>=0; h--) {
	 for(int g=nd-1; g>=0; g--) {
		double xs = -R1;		  double xe = -R1/4+(W1/2)/4;
		double ys = -R2 + (2*g-1)*W2/2;		double ye = -R2 + (2*g+3)*W2/2;
		double zs = -R3 + (2*h-1)*W3/2;		double ze = -R3 + (2*h+3)*W3/2;
		int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		double tx, ty, tz;
		tx = -R1/2;		  ty = (-R2 + (g+0.5)*W2)/2;		  tz = (-R3 + (h+0.5)*W3)/2;
		fx[wcnt] = tx;		  fy[wcnt] = ty;		  fz[wcnt] = tz;
		nx[wcnt] = xn;		  ny[wcnt] = yn;		  nz[wcnt] = zn;
		wcnt++;
	 }
  }

  //face 4: -y,-z,-x
  for(int f=nd-1; f>=0; f--) {
	 for(int h=nd-1; h>=0; h--) {
		double ys = -R2;		  double ye = -R2/4+(W2/2)/4;
		double zs = -R3 + (2*h-1)*W3/2;		  double ze = -R3 + (2*h+3)*W3/2;
		double xs = -R1 + (2*f-1)*W1/2;		  double xe = -R1 + (2*f+3)*W1/2;
		int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		double tx, ty, tz;
		ty = -R2/2;		  tz = (-R3 + (h+0.5)*W3)/2;		  tx = (-R1 + (f+0.5)*W1)/2;
		fx[wcnt] = tx;		  fy[wcnt] = ty;		  fz[wcnt] = tz;
		nx[wcnt] = xn;		  ny[wcnt] = yn;		  nz[wcnt] = zn;
		wcnt++;
	 }
  }
  
  //face 5: -z,-x,-y
  for(int g=nd-1; g>=0; g--) {
	 for(int f=nd-1; f>=0; f--) {
		double zs = -R3;		  double ze = -R3/4+(W3/2)/4;
		double xs = -R1 + (2*f-1)*W1/2;		double xe = -R1 + (2*f+3)*W1/2;
		double ys = -R2 + (2*g-1)*W2/2;		double ye = -R2 + (2*g+3)*W2/2;
		int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		double tx, ty, tz;
		tz = -R3/2;		  tx = (-R1 + (f+0.5)*W1)/2;		  ty = (-R2 + (g+0.5)*W2)/2;
		fx[wcnt] = tx;		  fy[wcnt] = ty;		  fz[wcnt] = tz;
		nx[wcnt] = xn;		  ny[wcnt] = yn;		  nz[wcnt] = zn;
		wcnt++;
	 }
  }
  
  return 0;
}

//----------------------------------------
int fdct3d_param(int N1, int N2, int N3, int nbscales, int nbdstz_coarse, int ac,
				  vector< vector<double> >& fxs, vector< vector<double> >& fys, vector< vector<double> >& fzs,
				  vector< vector<int   > >& nxs, vector< vector<int   > >& nys, vector< vector<int   > >& nzs)
{
  fxs.resize(nbscales);  fys.resize(nbscales);  fzs.resize(nbscales);
  nxs.resize(nbscales);  nys.resize(nbscales);  nzs.resize(nbscales);
  
  int L = nbscales;
  if(ac==1) {
	 {
		int s = 0;
		double L1 = 4.0*N1/3.0 / pow2(L-1-s);	 double L2 = 4.0*N2/3.0 / pow2(L-1-s);	 double L3 = 4.0*N3/3.0 / pow2(L-1-s);
		fdct3d_param_center(L1,L2,L3,s, fxs,fys,fzs, nxs,nys,nzs);
	 }
	 for(int s=1; s<L; s++) {
		double L1 = 4.0*N1/3.0 / pow2(L-1-s);	 double L2 = 4.0*N2/3.0 / pow2(L-1-s);	 double L3 = 4.0*N3/3.0 / pow2(L-1-s);
		int nd = nbdstz_coarse * pow2(s/2);
		fdct3d_param_angles(L1,L2,L3,s, nd, fxs,fys,fzs, nxs,nys,nzs);
	 }
  } else {
	 {
		int s = 0;
		double L1 = 4.0*N1/3.0 / pow2(L-1-s);	 double L2 = 4.0*N2/3.0 / pow2(L-1-s);	 double L3 = 4.0*N3/3.0 / pow2(L-1-s);
		fdct3d_param_center(L1,L2,L3,s, fxs,fys,fzs, nxs,nys,nzs);
	 }
	 for(int s=1; s<L-1; s++) {
		double L1 = 4.0*N1/3.0 / pow2(L-1-s);	 double L2 = 4.0*N2/3.0 / pow2(L-1-s);	 double L3 = 4.0*N3/3.0 / pow2(L-1-s);
		int nd = nbdstz_coarse * pow2(s/2);
		fdct3d_param_angles(L1,L2,L3,s, nd, fxs,fys,fzs, nxs,nys,nzs);
	 }
	 {
		int s = L-1;
		double L1 = 4.0*N1/3.0 / pow2(L-1-s);	 double L2 = 4.0*N2/3.0 / pow2(L-1-s);	 double L3 = 4.0*N3/3.0 / pow2(L-1-s);
		fdct3d_param_wavelet(L1,L2,L3,s, N1,N2,N3, fxs,fys,fzs, nxs,nys,nzs);
	 }
  }
  
  return 0;
}
