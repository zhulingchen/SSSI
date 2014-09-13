/*
  Copyright (C) 2004 Caltech
  Written by Lexing Ying
*/

#include "fdct_usfft.hpp"
#include "fdct_usfft_inline.hpp"

using namespace std;
FDCT_USFFT_NS_BEGIN_NAMESPACE

//------------------------------------
extern int fdct_usfft_adjfftL(int N1, int N2, CpxOffMat& T, CpxNumMat& x);
extern int fdct_usfft_adjsepscale(int N1, int N2, int nbscales, int ac, vector<CpxOffMat>& Xhghs, CpxOffMat& O);
extern int fdct_usfft_adjsepangle(double XL1, double XL2, int nbangles, vector<CpxOffMat>& msc, CpxOffMat& Xhgh);
extern int fdct_usfft_adjifftS(vector<CpxNumMat>& csc, vector<CpxOffMat>& msc);
extern int fdct_usfft_adjwavelet(vector<CpxNumMat>& csc, CpxOffMat& Xhgh);

extern int fdct_usfft_sepangle(double XL1, double XL2, int nbangles, CpxOffMat& Xhgh, vector<CpxOffMat>& msc);

int fdct_usfft_invsepangle(double XL1, double XL2, int nbangles, vector<CpxOffMat>& msc, CpxOffMat& Xhgh);
int fdct_usfft_makeToeplitz(double XL1, double XL2, int nbangles, vector< vector<CpxNumVec> >& Lambda,
			    map<int, fftw_plan>& f1map, map<int, fftw_plan>& b1map);
int fdct_usfft_multToeplitz(double XL1, double XL2, int nbangles, vector< vector<CpxNumVec> >& Lambda, CpxOffMat& d, CpxOffMat& p,
			    map<int, fftw_plan>& f1map, map<int, fftw_plan>& b1map);

//------------------------------------
int ifdct_usfft(int N1, int N2, int nbscales, int nbangles_coarse, int ac, vector< vector<CpxNumMat> >& c, CpxNumMat& x)
{
  assert(nbscales==c.size() && nbangles_coarse==c[1].size());
  //1.
  vector<CpxOffMat> Xhghs;  Xhghs.resize(nbscales);
  if(ac==1) {
    vector<int> nbangles(nbscales);
    nbangles[0] = 1;
    for(int sc=1; sc<nbscales; sc++)	 nbangles[sc] = nbangles_coarse * pow2( int(ceil(double(sc-1)/2)) );
    //finest+mid levels
    double XL1 = 4.0*N1/3.0;  double XL2 = 4.0*N2/3.0; //range
    for(int sc=nbscales-1; sc>0; sc--) {
      vector<CpxOffMat> msc(nbangles[sc]);
      fdct_usfft_adjifftS(c[sc], msc);
      fdct_usfft_invsepangle(XL1, XL2, nbangles[sc], msc, Xhghs[sc]);
      XL1 = XL1/2;	 XL2 = XL2/2;
    }
    fdct_usfft_adjwavelet(c[0], Xhghs[0]);
  } else {
    vector<int> nbangles(nbscales);
    nbangles[0] = 1;
    for(int sc=1; sc<nbscales-1; sc++)		nbangles[sc] = nbangles_coarse * pow2( int(ceil(double(sc-1)/2)) );
    nbangles[nbscales-1] = 1;
    //finest level
    fdct_usfft_adjwavelet(c[nbscales-1], Xhghs[nbscales-1]);
    //mid levels
    double XL1 = 2.0*N1/3.0;	 double XL2 = 2.0*N2/3.0;
    for(int sc=nbscales-2; sc>0; sc--) {
      vector<CpxOffMat> msc(nbangles[sc]);
      fdct_usfft_adjifftS(c[sc], msc);
      fdct_usfft_invsepangle(XL1, XL2, nbangles[sc], msc, Xhghs[sc]);
      XL1 = XL1/2;	 XL2 = XL2/2;
    }
    fdct_usfft_adjwavelet(c[0], Xhghs[0]);
  }
  //2.
  CpxOffMat O;
  fdct_usfft_adjsepscale(N1, N2, nbscales, ac, Xhghs, O);
  //3.
  fdct_usfft_adjfftL(N1, N2, O, x);
  return 0;
}

//--------------------------------------------
int fdct_usfft_invsepangle(double XL1, double XL2, int nbangles, vector<CpxOffMat>& msc, CpxOffMat& x)
{
  int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
  map<int, fftw_plan> f1map; //forward 1d map, IN_PLACE
  map<int, fftw_plan> b1map; //backward 1d map, IN_PLACE
  
  //Toeplitz preparation
  vector< vector<CpxNumVec> > Lambda;  fdct_usfft_makeToeplitz(XL1, XL2, nbangles, Lambda, f1map, b1map);

  //0. allocate space for x and set x = 0;
  x.resize(XS1, XS2);  clear(x);
  //1. call adjoint first, put into b
  CpxOffMat b(XS1, XS2);
  fdct_usfft_adjsepangle(XL1, XL2, nbangles, msc, b);
  //2. solve for MX = B, where M = (A^t A)
  
  CpxOffMat r(b);
  CpxOffMat d(r);
  CpxOffMat q(XS1, XS2); //empty
  
  cpx* rdata = r.data();  //cpx* xdata = x.data();  cpx* ddata = d.data();  cpx* qdata = q.data();
  
  int m = x.m();  int n = x.n();
  double v0=0;  for(int i=0; i<m*n; i++)	 v0 += norm(rdata[i]);
  double v1=v0;
  
  int it = 0;
  while(v1>1e-8 && it<25) {
    //q = Ad
    fdct_usfft_multToeplitz(XL1, XL2, nbangles, Lambda, d, q, f1map, b1map);
    
    cpx* rdata = r.data();	 cpx* xdata = x.data();	 cpx* ddata = d.data();	 cpx* qdata = q.data();
    //tmp = d' q = d' A d
    cpx tz(0,0);	 for(int i=0; i<m*n; i++)		tz += conj(ddata[i]) * qdata[i];
    double tmp=real(tz);
    
    //alpha = v1/tmp;
    double alpha = v1/tmp;
    //x = x + alpha d
    for(int i=0; i<m*n; i++)		xdata[i] += alpha * ddata[i];
    //r = r - alpha q
    for(int i=0; i<m*n; i++)		rdata[i] -= alpha * qdata[i];
    //reset v0 and v1
    v0 = v1;
    v1 = 0;	 for(int i=0; i<m*n; i++)		v1 += norm(rdata[i]);
    //beta = v1/v0;
    double beta = v1/v0;
    //d = r + beta d
    for(int i=0; i<m*n; i++)		ddata[i] = rdata[i] + beta * ddata[i];
    it++;
  }
  return 0;
}

//--------------------------------------------
int fdct_usfft_makeToeplitz(double XL1, double XL2, int nbangles, vector< vector<CpxNumVec> >& Lambda,
									 map<int, fftw_plan>& f1map, map<int, fftw_plan>& b1map)
{
  int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
  Lambda.resize(4); //4 quadrant
  int nd = nbangles / 4;
  int nbquadrants = 4;
  int wcnt = 0;

  double XL1b = XL1;  double XL2b = XL2;
  int qvec[] = {2,1,0,3};
  for(int qi=0; qi<nbquadrants; qi++) {
    int q = qvec[qi];
    //rotate
    fdct_usfft_rotate_forward(q, XL1b, XL2b, XL1, XL2);	 XL1 = abs(XL1);	 XL2 = abs(XL2);
    
    int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);	 //double XW1 = XL1/nd;
    double XW2 = XL2/nd;
    double xs = XR1/4; 		double xe = XR1;	 int xf = int(ceil(xs));
    
    int xn = int(ceil(xe-xs));	 int yn = 2*int(ceil(XW2))+1;
    if(xn%2==0) xn++;	 if(yn%2==0) yn++;	 //int xh = xn/2;
    int yh = yn/2;
    
    //allocate memory
    vector<CpxNumVec>& curLambda = Lambda[q];
    curLambda.resize( xn ); //xn columns
    
    for(int xcur=xf; xcur<xe; xcur++) { //for each line construct lambda
      int xid = xcur-xf;
      DblNumMat ysft(yn, nd);
      DblNumMat wsft(yn, nd);
      for(int w=0; w<nd; w++) {
	double ys = -XR2 + w*XW2-XW2/2;		  double ym = -XR2 + w*XW2+XW2/2;		  double ye = -XR2 + (w+1)*XW2+XW2/2;
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
      //construct vector
      for(int i=0; i<yn; i++)
	for(int j=0; j<nd; j++)
	  wsft(i,j) = wsft(i,j)*wsft(i,j); //squaring
      
      //get summation
      DblNumMat& off = ysft; //just reference
      DblNumMat& wgt = wsft; //just refernece
      int N = XS2;		int F = XF2;		assert(N%2==1);
      //fftw plan
      fftw_plan fp = NULL;
      map<int, fftw_plan>::iterator fit = f1map.find(N);
      if(fit!=f1map.end()) {		  fp = (*fit).second;
      } else {		  fp = fftw_create_plan(N, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);		  f1map[N] = fp;
      }
      CpxNumVec sum(N);
      CpxOffVec vvv(N);
      CpxNumVec feq(N);
      CpxNumVec tmp(N);
      for(int w=0; w<nd; w++) {
	for(int k=0; k<off.m(); k++) {
	  int rk = (k>F)? k-N : k;			 vvv(rk) = wgt(k,w);
	}
	fdct_usfft_ifftshift(vvv, tmp);
	fftw_one(fp, (fftw_complex*)tmp.data(), NULL);
	double scale = sqrt(double(N));	 for(int k=0; k<N; k++) {	 tmp(k) /= scale;  } //scaling
	for(int k=0; k<N; k++) {
	  double phase = -2*M_PI*k/double(N) * off(0,w);
	  feq(k) = tmp(k) * polar(1.0, phase);
	}
	for(int k=0; k<N; k++)
	  sum(k) += feq(k);
      } //all the wedges
      //divide it by sqrt(n)
      double scale = sqrt(double(N));		for(int k=0; k<N; k++)  		  sum(k) /= scale;
      //extend it
      CpxNumVec& tt = curLambda[xid];
      tt.resize(2*N-1);		clear(tt);
      for(int k=0; k<N; k++) {		  tt(k) = sum(k); }
      for(int k=1; k<N; k++) {		  tt(2*N-1-k) = conj(sum(k)); }
      //do ifft on it and store
      fftw_plan bp = NULL;
      map<int, fftw_plan>::iterator bit = b1map.find(2*N-1);
      if(bit!=b1map.end()) {		  bp = (*bit).second;
      } else {		  bp = fftw_create_plan(2*N-1, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);		  b1map[2*N-1] = bp;
      }
      fftw_one(bp, (fftw_complex*)tt.data(), NULL);
      scale = sqrt(double(2*N-1));		for(int k=0; k<2*N-1; k++)		  tt(k) /= scale; //scaling
    } //each line
    wcnt += nd;
  } //quadrant
  assert(wcnt==nbangles);
  return 0;
}

int fdct_usfft_multToeplitz(double XL1, double XL2, int nbangles, vector< vector<CpxNumVec> >& Lambda, CpxOffMat& d, CpxOffMat& z,
									 map<int, fftw_plan>& f1map, map<int, fftw_plan>& b1map)
{
  int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
  int nd = nbangles / 4;
  int nbquadrants = 4;
  int wcnt = 0;
  
  //clear space
  clear(z);
  CpxOffMat db(d);
  CpxOffMat zb(z);

  double XL1b = XL1;  double XL2b = XL2;
  int qvec[] = {2,1,0,3};
  for(int qi=0; qi<nbquadrants; qi++) {
	 int q = qvec[qi];
	 //rotate
	 fdct_usfft_rotate_forward(q, XL1b, XL2b, XL1, XL2);	 XL1 = abs(XL1);	 XL2 = abs(XL2);
	 fdct_usfft_rotate_forward(q, db, d);	 fdct_usfft_rotate_forward(q, zb, z);
	 
	 int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);	 //double XW1 = XL1/nd;
	 double XW2 = XL2/nd;
	 double xs = XR1/4; 		double xe = XR1;	 int xf = int(ceil(xs));
	 
	 int xn = int(ceil(xe-xs));	 int yn = 2*int(ceil(XW2))+1;
	 if(xn%2==0) xn++;	 if(yn%2==0) yn++;
	 //int xh = xn/2;	 int yh = yn/2;
	 
	 vector<CpxNumVec>& curLambda = Lambda[q];
	 
	 //for each line
	 for(int xcur=xf; xcur<xe; xcur++) { //for each line construct lambda
		int xid = xcur-xf;
		int N = XS2;		int F = XF2;		assert(N%2==1);
		
		map<int,fftw_plan>::iterator fit;
		fftw_plan sfp = NULL;
		fit=f1map.find(N);
		if(fit!=f1map.end()) {		  sfp = (*fit).second;
		} else {		  sfp = fftw_create_plan(N, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);		  f1map[N] = sfp;
		}
		fftw_plan lfp = NULL;
		fit=f1map.find(2*N-1);
		if(fit!=f1map.end()) {		  lfp = (*fit).second;
		} else {		  lfp = fftw_create_plan(2*N-1, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);		  f1map[2*N-1] = lfp;
		}
		map<int, fftw_plan>::iterator bit;
		fftw_plan sbp = NULL;
		bit=b1map.find(N);
		if(bit!=b1map.end()) {		  sbp = (*bit).second;
		} else {		  sbp = fftw_create_plan(N, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);		  b1map[N] = sbp;
		}
		fftw_plan lbp = NULL;
		bit=b1map.find(2*N-1);
		if(bit!=b1map.end()) {		  lbp = (*bit).second;
		} else {		  lbp = fftw_create_plan(2*N-1, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);		  b1map[2*N-1] = lbp;
		}
		
		CpxOffVec val(XS2);
		CpxNumVec tmp(N);
		CpxOffVec feq(N);
		CpxNumVec feqext(2*N-1);
		CpxNumVec newext(2*N-1);
		//1. get value and padding 0
		for(int i=-XF2; i<-XF2+XS2; i++) {		  val(i) = d(xcur,i); }
		//2. fft to get feq
		fdct_usfft_ifftshift(val, tmp); //shift
		fftw_one(sfp, (fftw_complex*)tmp.data(), NULL);
		double scale = sqrt(double(N));  for(int k=0; k<N; k++) {	 tmp(k) /= scale;	 } //scaling
		fdct_usfft_fftshift(tmp, feq); //tmp2 contains the right frequency
		//3. times with lambda
		CpxNumVec& tt = curLambda[xid];
		//   extend feq to feqext
		clear(feqext);
		for(int i=0; i<N; i++) {		  feqext(i) = feq(i-F); }
		//   large ifft
		fftw_one(lbp, (fftw_complex*)feqext.data(), NULL);
		scale = sqrt(double(2*N-1));		for(int i=0; i<2*N-1; i++) {		  feqext(i) /= scale; }
		//   multiplication
		for(int i=0; i<2*N-1; i++) {
		  newext(i) = feqext(i) * tt(i) * scale;
		}
		//   large fft
		fftw_one(lfp, (fftw_complex*)newext.data(), NULL);
		scale = sqrt(double(2*N-1));		for(int i=0; i<2*N-1; i++) {		  newext(i) /= scale; }
		//   cut feq
		for(int i=0; i<N; i++) {		  feq(i-F) = newext(i); }
		//4. fft back to get value
		fdct_usfft_ifftshift(feq, tmp); //shift
		fftw_one(sbp, (fftw_complex*)tmp.data(), NULL);
		scale = sqrt(double(N));		for(int k=0; k<N; k++) {		tmp(k) /= scale; }
		fdct_usfft_fftshift(tmp, val);
		//5. add back to z
		for(int i=-XF2; i<-XF2+XS2; i++) {		  z(xcur,i) += val(i); }
	 } //end of each line
	 wcnt += nd;
	 //rotate matrix
	 fdct_usfft_rotate_backward(q, z, zb);
  }
  
  //RESTORE THE VALUE
  z = zb;
  d = db;
  assert(wcnt==nbangles);
  return 0;
}

FDCT_USFFT_NS_END_NAMESPACE
