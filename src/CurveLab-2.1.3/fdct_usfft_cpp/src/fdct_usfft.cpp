/*
  Copyright (C) 2004 Caltech
  Written by Lexing Ying
*/

#include "fdct_usfft.hpp"
#include "fdct_usfft_inline.hpp"

FDCT_USFFT_NS_BEGIN_NAMESPACE

//------------------------------------
int fdct_usfft_fftL(int N1, int N2, CpxNumMat& x, CpxOffMat& O);
int fdct_usfft_sepscale(int N1, int N2, int nbscales, int ac, CpxOffMat& O, vector<CpxOffMat>& Xhghs);
int fdct_usfft_sepangle(double XL1, double XL2, int nbangles, CpxOffMat& Xhgh, vector<CpxOffMat>& msc);
int fdct_usfft_ifftS(vector<CpxOffMat>& msc, vector<CpxNumMat>& csc);
int fdct_usfft_wavelet(CpxOffMat& Xhgh, vector<CpxNumMat>& csc);

int fdct_usfft_1dinterp(DblNumMat& off, DblNumMat& wgt, CpxOffVec& val, CpxNumMat& res,
								map<int, fftw_plan>& f1map, map<int, fftw_plan>& b1map);

//------------------------------------
int fdct_usfft( int N1, int N2, int nbscales, int nbangles_coarse, int ac, CpxNumMat& x, vector< vector<CpxNumMat> >& c)
{
  assert(N1==x.m() && N2==x.n());
  int F1 = N1/2;  int F2 = N2/2;
  
  //1. fftL
  CpxOffMat O(N1, N2);  fdct_usfft_fftL(N1, N2, x, O);
  
  //2. seperate scale
  vector<CpxOffMat> Xhghs(nbscales);
  fdct_usfft_sepscale(N1, N2, nbscales, ac, O, Xhghs);
  
  //3. work on each scale
  if(ac==1) {
	 //nbangles
	 vector<int> nbangles(nbscales);
	 nbangles[0] = 1;
	 for(int sc=1; sc<nbscales; sc++)		nbangles[sc] = nbangles_coarse * pow2( int(ceil(double(sc-1)/2)) );
	 //c
	 c.resize(nbscales);
	 for(int sc=0; sc<nbscales; sc++)		c[sc].resize( nbangles[sc] );
	 //finest+mid levels
	 double XL1 = 4.0*N1/3.0;  double XL2 = 4.0*N2/3.0; //range
	 for(int sc=nbscales-1; sc>0; sc--) {
		vector<CpxOffMat> msc(nbangles[sc]);
		fdct_usfft_sepangle(XL1, XL2, nbangles[sc], Xhghs[sc], msc);
		fdct_usfft_ifftS(msc, c[sc]);
		XL1 = XL1/2;	 XL2 = XL2/2;
	 }
	 //coarest level
	 fdct_usfft_wavelet(Xhghs[0], c[0]);
  } else {
	 //nbangles
	 vector<int> nbangles(nbscales);
	 nbangles[0] = 1;
	 for(int sc=1; sc<nbscales-1; sc++)		nbangles[sc] = nbangles_coarse * pow2( int(ceil(double(sc-1)/2)) );
	 nbangles[nbscales-1] = 1;
	 //c
	 c.resize(nbscales);
	 for(int sc=0; sc<nbscales; sc++)		c[sc].resize( nbangles[sc] );
	 //finest level
	 fdct_usfft_wavelet(Xhghs[nbscales-1], c[nbscales-1]);
	 //mid levels
	 double XL1 = 2.0*N1/3.0;  double XL2 = 2.0*N2/3.0; //range
	 for(int sc=nbscales-2; sc>0; sc--) {
		vector<CpxOffMat> msc(nbangles[sc]);
		fdct_usfft_sepangle(XL1, XL2, nbangles[sc], Xhghs[sc], msc);
		fdct_usfft_ifftS(msc, c[sc]);
		XL1 = XL1/2;	 XL2 = XL2/2;
	 }
	 //coarest level
	 fdct_usfft_wavelet(Xhghs[0], c[0]);
  }
  return 0;
}

//--------------------
int fdct_usfft_fftL(int N1, int N2, CpxNumMat& x, CpxOffMat& O)
{
  //move zero freq at front
  int F1 = N1/2;  int F2 = N2/2;
  //1. ifftshift
  CpxNumMat T(x);
  //2. do fft and scale
  fftwnd_plan p = fftw2d_create_plan(N2, N1, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fftwnd_one(p, (fftw_complex*)T.data(), NULL);
  fftwnd_destroy_plan(p);
  double sqrtprod = sqrt(double(N1*N2));  for(int i=0; i<N1; i++)	 for(int j=0; j<N2; j++)		T(i,j) /= sqrtprod;
  //3. fftshift
  O.resize(N1, N2);
  fdct_usfft_fftshift(T, O);
  return 0;
}

//--------------------
int fdct_usfft_sepscale(int N1, int N2, int nbscales, int ac, CpxOffMat& O, vector<CpxOffMat>& Xhghs)
{
  //-----------------------------------------------------
  //unfold if necessary
  CpxOffMat X;
  if(ac==1) {
	 //--------------------------
	 double XL1 = 4.0*N1/3.0;  double XL2 = 4.0*N2/3.0; //range
	 int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
	 IntOffVec t1(XS1);
	 for(int i=-XF1; i<-XF1+XS1; i++)		if(     i<-N1/2) t1(i) = i+int(N1);		else if(i>(N1-1)/2) t1(i) = i-int(N1);		else t1(i) = i;						 
	 IntOffVec t2(XS2);
	 for(int i=-XF2; i<-XF2+XS2; i++)		if(     i<-N2/2) t2(i) = i+int(N2);		else if(i>(N2-1)/2) t2(i) = i-int(N2);		else t2(i) = i;
	 X.resize(XS1, XS2);
	 for(int i=-XF1; i<-XF1+XS1; i++)
		for(int j=-XF2; j<-XF2+XS2; j++)
		  X(i,j) = O(t1(i), t2(j));
	 DblOffMat lowpass(XS1,XS2);
	 fdct_usfft_lowpasscompute(XL1, XL2, lowpass); //compute the low pass filter
	 for(int i=-XF1; i<-XF1+XS1; i++)
		for(int j=-XF2; j<-XF2+XS2; j++)
		  X(i,j) *= lowpass(i,j);
  } else {
	 //--------------------------
	 X = O;
  }
  //-----------------------------------------------------
  //seperate
  Xhghs.resize(nbscales);
  double XL1 = 4.0*N1/3.0;  double XL2 = 4.0*N2/3.0; //range
  //int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
  for(int sc=nbscales-1; sc>0; sc--) {
	 double XL1n = XL1/2;	 double XL2n = XL2/2;
	 int XS1n, XS2n;  int XF1n, XF2n;  double XR1n, XR2n;
	 fdct_usfft_rangecompute(XL1n, XL2n, XS1n, XS2n, XF1n, XF2n, XR1n, XR2n);
	 //get filters
	 DblOffMat lowpass(XS1n, XS2n);
	 fdct_usfft_lowpasscompute(XL1n, XL2n, lowpass);
	 DblOffMat hghpass(XS1n, XS2n);
	 for(int i=-XF1n; i<-XF1n+XS1n; i++)
		for(int j=-XF2n; j<-XF2n+XS2n; j++)
		  hghpass(i,j) = sqrt(1-lowpass(i,j)*lowpass(i,j));
	 //get Xhgh
	 CpxOffMat Xhgh(X);
	 for(int i=-XF1n; i<-XF1n+XS1n; i++)
		for(int j=-XF2n; j<-XF2n+XS2n; j++) {		  //Xhgh(i,j).re *= hghpass(i,j);		  Xhgh(i,j).im *= hghpass(i,j);
		  Xhgh(i,j) *= hghpass(i,j);
		}
	 CpxOffMat Xlow(XS1n, XS2n);
	 for(int i=-XF1n; i<-XF1n+XS1n; i++)
		for(int j=-XF2n; j<-XF2n+XS2n; j++) {		  //Xlow(i,j).re = X(i,j).re * lowpass(i,j);		  Xlow(i,j).im = X(i,j).im * lowpass(i,j);
		  Xlow(i,j) = X(i,j) * lowpass(i,j);
		}
	 //set into vector
	 Xhghs[sc] = Xhgh;
	 X = Xlow;
	 XL1 = XL1/2;	 XL2 = XL2/2;
  }
  Xhghs[0] = X;
  return 0;
}

//--------------------
int fdct_usfft_sepangle(double XL1, double XL2, int nbangles, CpxOffMat& Xhgh, vector<CpxOffMat>& msc)
{
  //int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
  map<int, fftw_plan> f1map; //forward 1d map, IN_PLACE
  map<int, fftw_plan> b1map; //backward 1d map, IN_PLACE
  
  //allocate msc
  msc.resize(nbangles);
  int nbquadrants = 4;
  int nd = nbangles / 4;
  int wcnt = 0;
  
  CpxOffMat Xhghb(Xhgh);
  double XL1b = XL1;  double XL2b = XL2;
  
  int qvec[] = {2,1,0,3};
  for(int qi=0; qi<nbquadrants; qi++) {
	 int q = qvec[qi];
	 //rotate data
  	 fdct_usfft_rotate_forward(q, XL1b, XL2b, XL1, XL2);	 XL1 = abs(XL1);	 XL2 = abs(XL2);
	 fdct_usfft_rotate_forward(q, Xhghb, Xhgh);
	 //  sample using USFFT for each line
	 int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);	 //double XW1 = XL1/nd;
	 double XW2 = XL2/nd;
	 double xs = XR1/4;		double xe = XR1;	 int xf = int(ceil(xs));
	 //xn,yn, make them odd
	 int xn = int(ceil(xe-xs));	 int yn = 2*int(ceil(XW2))+1;
	 if(xn%2==0) xn++;	 if(yn%2==0) yn++;
	 int xh = xn/2;	 int yh = yn/2;
	 //allocate temporary space for tsc
	 vector<CpxOffMat> tsc(nd);
	 for(int w=0; w<nd; w++)		tsc[w].resize(xn,yn);
	 //gather weighted samples for each wedges
	 for(int xcur=xf; xcur<xe; xcur++) { //for each line
		DblNumMat ysft(yn, nd); //store shifted positions
		DblNumMat wsft(yn, nd); //store shifted weights
		for(int w=0; w<nd; w++) {
		  double ys = -XR2 + (w-0.5)*XW2;		  double ym = -XR2 + (w+0.5)*XW2;		  double ye = -XR2 + (w+1.5)*XW2;
		  double s0 = ys/XR1;		  double s2 = ym/XR1;		  double s4 = ye/XR1;
		  for(int yid=0; yid<yn; yid++) {
			 double tmp = (yid-yh) + s2*xcur; //shifting operator, 
			 double pou;
			 if(tmp<s2*xcur) { //below
				double l,r; fdct_usfft_window( (tmp/xcur-s0)/(s2-s0), l, r);				  pou = l;
			 } else {
				double l,r; fdct_usfft_window( (tmp/xcur-s2)/(s4-s2), l, r);				  pou = r;
			 }
			 ysft(yid,w) = tmp;
			 wsft(yid,w) = pou;
		  }
		}
		//interpolation with weights (usfft or usfft_simple)
		CpxOffVec val(XS2); 		////CpxOffVec val(XS2+1, -XF2-1);  //value for the current column, the first element is zero
		for(int i=-XF2; i<-XF2+XS2; i++)		  val(i) = Xhgh(xcur,i);
		CpxNumMat res(yn, nd);
		fdct_usfft_1dinterp(ysft, wsft, val, res, f1map, b1map); //with weights 
		
		int tmpx = xcur%xn;		if(tmpx<-xh) tmpx+=xn;				  if(tmpx>=-xh+xn) tmpx-=xn;
		for(int w=0; w<nd; w++)
		  for(int yid=0; yid<yn; yid++)
			 tsc[w](tmpx,yid-yh) = res(yid,w);
	 }//each line
	 
	 //rotate data back into msc
	 for(int w=nd-1; w>=0; w--) {
		fdct_usfft_rotate_backward(q, tsc[w], msc[wcnt]);
		wcnt++;
	 }
  } //for loop for quadrant
  Xhgh = Xhghb;
  XL1 = XL1b;  XL2 = XL2b;
  assert(wcnt==nbangles);
  
  for(map<int,fftw_plan>::iterator mit=f1map.begin(); mit!=f1map.end(); mit++) {
	 fftw_plan p = (*mit).second;
	 fftw_destroy_plan(p);
  }
  for(map<int,fftw_plan>::iterator mit=b1map.begin(); mit!=b1map.end(); mit++) {
	 fftw_plan p = (*mit).second;
	 fftw_destroy_plan(p);
  }
  return 0;
}

//--------------------
int fdct_usfft_ifftS(vector<CpxOffMat>& msc, vector<CpxNumMat>& csc)
{
  typedef pair<int,int> intpair;
  map<intpair, fftwnd_plan> planmap;
  //do work
  csc.resize(msc.size());
  for(int w=0; w<msc.size(); w++) {
	 //allocate space
	 int xn = msc[w].m();	 int yn = msc[w].n();
	 int xh = xn/2;	 int yh = yn/2;
	 //shift
	 CpxNumMat tpdata(xn,yn);	 fdct_usfft_ifftshift(msc[w],tpdata);
	 //fft
	 map<intpair,fftwnd_plan>::iterator mit=planmap.find( intpair(xn,yn) );
	 fftwnd_plan p = NULL;
	 if(mit!=planmap.end()) {
		p = (*mit).second;
	 } else {
		p = fftw2d_create_plan(yn, xn, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
		planmap[ intpair(xn,yn) ] = p;
	 }
	 fftwnd_one(p, (fftw_complex*)tpdata.data(), NULL);
	 double sqrtprod = sqrt(double(xn*yn));
	 for(int i=0; i<xn; i++)		for(int j=0; j<yn; j++)		  tpdata(i,j) /= sqrtprod;
	 //shift
	 csc[w] = tpdata;	 //csc[w].resize(xn,yn,-xh,-yh);	 fdct_usfft_fftshift(xn,yn,xh,yh,tpdata,csc[w]);
  }
  //delete planners
  for(map<intpair,fftwnd_plan>::iterator mit=planmap.begin(); mit!=planmap.end(); mit++) {
	 fftwnd_plan p = (*mit).second;
	 fftwnd_destroy_plan(p);
  }
  return 0;
}

//--------------------
int fdct_usfft_wavelet(CpxOffMat& Xhgh, vector<CpxNumMat>& csc)
{
  int N1 = Xhgh.m();  int N2 = Xhgh.n();
  int F1 = -Xhgh.s();  int F2 = -Xhgh.t();
  CpxNumMat T(N1, N2);
  fdct_usfft_ifftshift(Xhgh, T);
  fftwnd_plan p = fftw2d_create_plan(N2, N1, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fftwnd_one(p, (fftw_complex*)T.data(), NULL);
  fftwnd_destroy_plan(p);
  double sqrtprod = sqrt(double(N1*N2));  for(int i=0; i<N1; i++)	 for(int j=0; j<N2; j++)		T(i,j) /= sqrtprod;
  csc[0] = T;   //csc[0].resize(N1, N2, -F1, -F2);  //fdct_usfft_fftshift(N1, N2, F1, F2, T, csc[0]);
  return 0;
}

//--------------------
int fdct_usfft_1dinterp(DblNumMat& off, DblNumMat& wgt, CpxOffVec& val, CpxNumMat& res,
								map<int, fftw_plan>& f1map, map<int, fftw_plan>& b1map)
{
  if(off.n()<=64) { //SIMPLE INTERPOLATION
	 //--------------------------------------------
	 int N = val.m();
	 int F = -val.s();
	 fftw_plan fp = NULL;
	 map<int, fftw_plan>::iterator fit = f1map.find(N);
	 if(fit!=f1map.end()) {	 fp = (*fit).second;
	 } else {	 fp = fftw_create_plan(N, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);	 f1map[N] = fp;
	 }
	 fftw_plan bp = NULL;
	 map<int, fftw_plan>::iterator bit = b1map.find(N);
	 if(bit!=b1map.end()) {	 bp = (*bit).second;
	 } else {	 bp = fftw_create_plan(N, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);	 b1map[N] = bp;
	 }
	 //CpxOffVec val //right value
	 CpxOffVec feq(N); //right freq
	 CpxOffVec vvv(N); //new value
	 CpxOffVec fff(N); //new freq
	 CpxNumVec tmp(N);
	 //1. fft_mid val
	 fdct_usfft_ifftshift(val, tmp);
	 fftw_one(fp, (fftw_complex*)tmp.data(), NULL); //fft
	 double scale = sqrt(double(N));  for(int k=0; k<N; k++)		tmp(k) /= scale; //scaling
	 fdct_usfft_fftshift(tmp, feq); //frequency
	 //2. for each wedge, multiply, ifft and put to res 
	 int nbwedges = off.n();
	 for(int w=0; w<nbwedges; w++) {
		for(int k=-F; k<-F+N; k++) {
		  double phase = 2.0*M_PI*double(k)/double(N) * off(0,w);
		  fff(k) = feq(k) * polar(1.0, phase);
		} //new frequency
		fdct_usfft_ifftshift(fff, tmp);
		fftw_one(bp, (fftw_complex*)tmp.data(), NULL);
		double scale = sqrt(double(N));	 for(int k=0; k<N; k++)		  tmp(k) /= scale; //scaling
		fdct_usfft_fftshift(tmp, vvv); //new value
		for(int k=0; k<off.m(); k++) {
		  int rk = (k>F) ? k-N : k;
		  res(k,w) = vvv(rk) * wgt(k,w);
		}
	 }
  } else { //USFFT INTERPOLATION
	 //--------------------------------------------
	 int L = 4;
	 int D = 16;
	 int N = val.m();
	 int F = -val.s();
	 fftw_plan fp = NULL;
	 map<int, fftw_plan>::iterator fit=f1map.find(N);
	 if(fit!=f1map.end()) {		fp = (*fit).second;
	 } else {	 fp = fftw_create_plan(N, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);	 f1map[N] = fp;
	 }
	 fftw_plan bp = NULL;
	 map<int, fftw_plan>::iterator bit = b1map.find(D*N);
	 if(bit!=b1map.end()) {	 bp = (*bit).second;
	 } else {	 bp = fftw_create_plan(D*N, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);	 b1map[D*N] = bp;
	 }
	 CpxNumVec tmp(N);
	 CpxNumVec exttmp(D*N);
	 CpxOffVec feq(N);
	 vector<CpxOffVec> extdat(L);
	 //1. get freq
	 fdct_usfft_ifftshift(val, tmp);
	 fftw_one(fp, (fftw_complex*)tmp.data(), NULL); //here tmp are the real frequencies
	 double scale = double(N);	 for(int k=0; k<N; k++) {	 tmp(k) /= scale;	 } //scaling
	 fdct_usfft_fftshift(tmp, feq); //frequency
	 //2. extend frequency and get frequencies of derivatives
	 extdat[0].resize(D*N);
	 for(int k=-F; k<-F+N; k++) {		extdat[0](k) = feq(k); }
	 for(int l=1; l<L; l++) {
		extdat[l].resize(D*N);
		for(int k=-F; k<-F+N; k++) {
		  extdat[l](k) = extdat[l-1](k) * cpx(0,k);
		}
	 }
	 //3. get derivatives
	 for(int l=0; l<L; l++) {
		fdct_usfft_ifftshift(extdat[l], exttmp);
		fftw_one(bp, (fftw_complex*)exttmp.data(), NULL);
		fdct_usfft_fftshift(exttmp, extdat[l]);
	 }	 //extdat[l] contains the lth derivative
	 //4. do tayler for each point
	 double step = 2.0* M_PI / double(D*N);
	 for(int j=0; j<off.n(); j++)
		for(int i=0; i<off.m(); i++) {
		  //1. find the grid point on the left, 
		  double cof = off(i,j) / double(N) * 2.0*M_PI;
		  if(cof<-M_PI)			 cof+= 2.0*M_PI;		  else if(cof>=M_PI)			 cof-= 2.0*M_PI; //collapse into [-pi,pi)
		  int ind = int(floor(cof/step));		  ind = min(max(ind, -D*N/2), D*N/2-1);
		  double dta = cof-ind*step;
		  double pow = 1;
		  cpx sum(0.0,0.0);//
		  for(int l=0; l<L; l++) {
			 sum += extdat[l](ind) * pow;
			 pow = pow * dta/(l+1);
		  }
		  res(i,j) = sum * wgt(i,j);
		}
  }
  return 0;
}

FDCT_USFFT_NS_END_NAMESPACE
