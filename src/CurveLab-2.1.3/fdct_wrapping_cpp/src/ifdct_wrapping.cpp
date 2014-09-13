/*
   Copyright (C) 2004 Caltech
   Written by Lexing Ying
*/

#include "fdct_wrapping.hpp"
#include "fdct_wrapping_inline.hpp"

FDCT_WRAPPING_NS_BEGIN_NAMESPACE

int fdct_wrapping_invsepangle(double XL1, double XL2, int nbangle, vector<CpxNumMat>& csc, CpxOffMat& Xhgh);
int fdct_wrapping_invwavelet(vector<CpxNumMat>& csc, CpxOffMat& Xhgh);

//-------------------------------------------------------------------------------
int ifdct_wrapping(int N1, int N2, int nbscales, int nbangles_coarse, int allcurvelets, vector< vector<CpxNumMat> >& c, CpxNumMat& x)
{
  assert(nbscales==c.size() && nbangles_coarse==c[1].size());
  
  int F1 = N1/2;  int F2 = N2/2;
  //-------------------------------------------angles to Xhgh
  vector<int> nbangles(nbscales);
  vector<CpxOffMat> Xhghs;  Xhghs.resize(nbscales);
  
  if(allcurvelets==1) {
	 //----
	 nbangles[0] = 1;
	 for(int sc=1; sc<nbscales; sc++)		nbangles[sc] = nbangles_coarse * pow2( int(ceil(double(sc-1)/2)) );
	 
	 double XL1 = 4.0*N1/3.0;  double XL2 = 4.0*N2/3.0;
	 for(int sc=nbscales-1; sc>0; sc--) {
		fdct_wrapping_invsepangle(XL1, XL2, nbangles[sc], c[sc], Xhghs[sc]);
		XL1 /= 2;		XL2 /= 2;
	 }
	 fdct_wrapping_invwavelet(c[0], Xhghs[0]);
  } else {
	 //----
	 nbangles[0] = 1;
	 for(int sc=1; sc<nbscales-1; sc++)		nbangles[sc] = nbangles_coarse * pow2( int(ceil(double(sc-1)/2)) );
	 nbangles[nbscales-1] = 1;
	 
	 fdct_wrapping_invwavelet(c[nbscales-1], Xhghs[nbscales-1]);
	 double XL1 = 2.0*N1/3.0;	 double XL2 = 2.0*N2/3.0;
	 for(int sc=nbscales-2; sc>0; sc--) {
		fdct_wrapping_invsepangle(XL1, XL2, nbangles[sc], c[sc], Xhghs[sc]);
		XL1 /= 2;		XL2 /= 2;
	 }
	 fdct_wrapping_invwavelet(c[0], Xhghs[0]);
  }
  
  //-------------------------------------------xhghs to O
  //combine
  CpxOffMat X;
  if(allcurvelets==1) {
	 double XL1 = 4.0*N1/3.0;  double XL2 = 4.0*N2/3.0; //range
	 int XS1, XS2;  int XF1, XF2;  double XR1, XR2;	 fdct_wrapping_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
	 X.resize(XS1, XS2);
  } else {
	 X.resize(N1, N2);
  }
  double XL1 = 4.0*N1/3.0;  double XL2 = 4.0*N2/3.0;
  int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_wrapping_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
  for(int sc=nbscales-1; sc>0; sc--) {
	 double XL1n = XL1/2;	 double XL2n = XL2/2;
	 int XS1n, XS2n;	 int XF1n, XF2n;	 double XR1n, XR2n;
	 fdct_wrapping_rangecompute(XL1n, XL2n, XS1n, XS2n, XF1n, XF2n, XR1n, XR2n);
	 
	 DblOffMat lowpass(XS1n, XS2n);
	 fdct_wrapping_lowpasscompute(XL1n, XL2n, lowpass);
	 DblOffMat hghpass(XS1n, XS2n);
	 for(int j=-XF2n; j<-XF2n+XS2n; j++)
		for(int i=-XF1n; i<-XF1n+XS1n; i++)
		hghpass(i,j) = sqrt(1-lowpass(i,j)*lowpass(i,j));
	 for(int j=-XF2n; j<-XF2n+XS2n; j++)
		for(int i=-XF1n; i<-XF1n+XS1n; i++)
		  Xhghs[sc](i,j) *= hghpass(i,j);
	 for(int j=-XF2n; j<-XF2n+XS2n; j++)
		for(int i=-XF1n; i<-XF1n+XS1n; i++)
		  Xhghs[sc-1](i,j) *= lowpass(i,j);
	 CpxOffMat& G = Xhghs[sc];
	 for(int j=G.t(); j<G.t()+G.n(); j++)
		for(int i=G.s(); i<G.s()+G.m(); i++)
		  X(i,j) += G(i,j);
	 XL1 = XL1/2;	 XL2 = XL2/2;
	 fdct_wrapping_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
  }
  for(int j=-XF2; j<-XF2+XS2; j++)
	 for(int i=-XF1; i<-XF1+XS1; i++)
		X(i,j) += Xhghs[0](i,j);
  // fold
  CpxOffMat O(N1, N2);
  if(allcurvelets==1) {
	 double XL1 = 4.0*N1/3.0;  double XL2 = 4.0*N2/3.0;
	 int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_wrapping_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
	 //times pou;
	 DblOffMat lowpass(XS1,XS2);
	 fdct_wrapping_lowpasscompute(XL1, XL2, lowpass);
	 for(int j=-XF2; j<-XF2+XS2; j++)
		for(int i=-XF1; i<-XF1+XS1; i++)
		  X(i,j) *= lowpass(i,j);
	 IntOffVec t1(XS1);
	 for(int i=-XF1; i<-XF1+XS1; i++)		if(     i<-N1/2) t1(i) = i+int(N1);		else if(i>(N1-1)/2) t1(i) = i-int(N1);		else t1(i) = i;
	 IntOffVec t2(XS2);
	 for(int i=-XF2; i<-XF2+XS2; i++)		if(     i<-N2/2) t2(i) = i+int(N2);		else if(i>(N2-1)/2) t2(i) = i-int(N2);		else t2(i) = i;
	 for(int j=-XF2; j<-XF2+XS2; j++)
		for(int i=-XF1; i<-XF1+XS1; i++)
		  O(t1(i), t2(j)) += X(i,j);
  } else {
	 O = X;
  }
  //------------------------------------------------------------
  CpxNumMat T(N1,N2);
  fdct_wrapping_ifftshift(O, T);
  fftwnd_plan p = fftw2d_create_plan(N2, N1, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fftwnd_one(p, (fftw_complex*)T.data(), NULL);
  fftwnd_destroy_plan(p);
  double sqrtprod = sqrt(double(N1*N2)); //scale
  for(int i=0; i<N1; i++)	 for(int j=0; j<N2; j++)	 T(i,j) /= sqrtprod;

  x = T;
  //x.resize(N1, N2);
  //fdct_wrapping_fftshift(T, x);
  
  return 0; 
}

//---------------------
int fdct_wrapping_invsepangle(double XL1, double XL2, int nbangle, vector<CpxNumMat>& csc, CpxOffMat& Xhgh)
{
  typedef pair<int,int> intpair;
  map<intpair, fftwnd_plan> planmap;
  
  int XS1, XS2;  int XF1, XF2;  double XR1, XR2;	 fdct_wrapping_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
  Xhgh.resize(XS1, XS2);
  
  int nbquadrants = 4;
  int nd = nbangle / 4;
  int wcnt = 0;
  
  //backup
  CpxOffMat Xhghb(Xhgh);
  double XL1b = XL1;  double XL2b = XL2;

  int qvec[] = {2,1,0,3};
  for(int qi=0; qi<nbquadrants; qi++) {
	 int q = qvec[qi];
	 //ROTATE data to its right position
	 fdct_wrapping_rotate_forward(q, XL1b, XL2b, XL1, XL2);	 XL1 = abs(XL1);	 XL2 = abs(XL2);
	 fdct_wrapping_rotate_forward(q, Xhghb, Xhgh);
	 //figure out XS, XF, XR
	 double XW1 = XL1/nd;	 double XW2 = XL2/nd;
	 int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_wrapping_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
	 for(int w=nd-1; w>=0; w--) {
		double xs = XR1/4 - (XW1/2)/4;
		double xe = XR1;
		double ys = -XR2 + (w-0.5)*XW2;
		double ye = -XR2 + (w+1.5)*XW2; //x range
		int xn = int(ceil(xe-xs));			 int yn = int(ceil(ye-ys));
		//MAKE THEM ODD
		if(xn%2==0) xn++;		if(yn%2==0) yn++;
		int xf = int(ceil(xs));			 //int yf = int(ceil(ys));
		//theta
		double thts, thtm, thte; //y direction
		if(w==0) {
		  thts = atan2(-1.0, 1.0-1.0/nd);
		  thtm = atan2(-1.0+1.0/nd, 1.0);
		  thte = atan2(-1.0+3.0/nd, 1.0);
		} else if(w==nd-1) {
		  thts = atan2(-1.0+(2.0*w-1.0)/nd, 1.0);
		  thtm = atan2(-1.0+(2.0*w+1.0)/nd, 1.0);
		  thte = atan2(1.0, 1.0-1.0/nd);
		} else {
		  thts = atan2(-1.0+(2.0*w-1.0)/nd, 1.0);
		  thtm = atan2(-1.0+(2.0*w+1.0)/nd, 1.0);
		  thte = atan2(-1.0+(2.0*w+3.0)/nd, 1.0);
		}
		int xh = xn/2;		int yh = yn/2; //half length
		CpxOffMat wpdata(xn,yn);

		{
		  //load
		  int xn = csc[wcnt].m();		  int yn = csc[wcnt].n();
		  CpxNumMat tpdata(csc[wcnt]);
		  //fft
		  fftwnd_plan p = NULL;
		  map<intpair, fftwnd_plan>::iterator mit=planmap.find( intpair(xn,yn) );
		  if(mit!=planmap.end()) {
			 p = (*mit).second;
		  } else {
			 p = fftw2d_create_plan(yn, xn, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
			 planmap[ intpair(xn, yn) ] = p;
		  }
		  fftwnd_one(p, (fftw_complex*)tpdata.data(), NULL);
		  double sqrtprod = sqrt(double(xn*yn));
		  for(int i=0; i<xn; i++)		  for(int j=0; j<yn; j++)			 tpdata(i,j) /= sqrtprod;
		  //fftshift
		  CpxOffMat rpdata;
		  fdct_wrapping_fftshift(tpdata,rpdata);
		  //rotate forward
		  fdct_wrapping_rotate_forward(q, rpdata, wpdata);
		}
		
		double R21 = XR2/XR1; //ratio
		for(int xcur=xf; xcur<xe; xcur++) { //for each layer
		  int yfm = (int)ceil( max(-XR2, R21*xcur*tan(thts)) );
		  int yto = (int)floor( min(XR2, R21*xcur*tan(thte)) );
		  for(int ycur=yfm; ycur<=yto; ycur++) {
			 int tmpx = xcur%xn;				  if(tmpx<-xh) tmpx+=xn;				  if(tmpx>=-xh+xn) tmpx-=xn;
			 int tmpy = ycur%yn;				  if(tmpy<-yh) tmpy+=yn;				  if(tmpy>=-yh+yn) tmpy-=yn;
			 //partition of unity
			 double thtcur = atan2(ycur/XR2, xcur/XR1);
			 double wtht;
			 if(thtcur<thtm) {
				double l,r; fdct_wrapping_window((thtcur-thts)/(thtm-thts), l, r);
				wtht = l;
			 } else {
				double l,r; fdct_wrapping_window((thtcur-thtm)/(thte-thtm), l, r);
				wtht = r;
			 }
			 double pou = wtht;
			 wpdata(tmpx,tmpy) *= pou;
			 Xhgh(xcur,ycur) += wpdata(tmpx,tmpy);
		  }
		}
		
		wcnt++;
	 }//w loop
	 fdct_wrapping_rotate_backward(q, Xhgh, Xhghb);
  } //q loop
  Xhgh = Xhghb;
  XL1 = XL1b;  XL2 = XL2b;

  assert(wcnt==nbangle);
  
  for(map<intpair, fftwnd_plan>::iterator mit=planmap.begin(); mit!=planmap.end(); mit++) {
	 fftwnd_plan p = (*mit).second;
	 fftwnd_destroy_plan(p);
  }
  return 0;
}

//---------------------
int fdct_wrapping_invwavelet(vector<CpxNumMat>& csc, CpxOffMat& Xhgh)
{
  assert(csc.size()==1);
  CpxNumMat& C = csc[0];
  int N1 = C.m();  int N2 = C.n();
  
  CpxNumMat T(C);  //CpxNumMat T(N1, N2);  fdct_wrapping_ifftshift(N1, N2, F1, F2, C, T);
  
  fftwnd_plan p = fftw2d_create_plan(N2, N1, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fftwnd_one(p, (fftw_complex*)T.data(), NULL);
  fftwnd_destroy_plan(p);
  double sqrtprod = sqrt(double(N1*N2));
  for(int j=0; j<N2; j++)
	 for(int i=0; i<N1; i++)
		T(i,j) /= sqrtprod;
  
  Xhgh.resize(N1, N2);
  fdct_wrapping_fftshift(T, Xhgh);
  return 0;
}

FDCT_WRAPPING_NS_END_NAMESPACE
