/*
  Copyright (C) 2004 Caltech
  Written by Lexing Ying
*/

#include "fdct_usfft.hpp"
#include "fdct_usfft_inline.hpp"

FDCT_USFFT_NS_BEGIN_NAMESPACE

//------------------------------------
int fdct_usfft_adjfftL(int N1, int N2, CpxOffMat& O, CpxNumMat& x);
int fdct_usfft_adjsepscale(int N1, int N2, int nbscales, int ac, vector<CpxOffMat>& Xhghs, CpxOffMat& O);
int fdct_usfft_adjsepangle(double XL1, double XL2, int nbangles, vector<CpxOffMat>& msc, CpxOffMat& Xhgh);
int fdct_usfft_adjifftS(vector<CpxNumMat>& csc, vector<CpxOffMat>& msc);
int fdct_usfft_adjwavelet(vector<CpxNumMat>& csc, CpxOffMat& Xhgh);

int fdct_usfft_adj1dinterp(DblNumMat& off, DblNumMat& wgt, CpxNumMat& res, CpxOffVec& val,
									map<int, fftw_plan>& f1map, map<int, fftw_plan>& b1map);
  
//------------------------------------
int afdct_usfft(int N1, int N2, int nbscales, int nbangles_coarse, int ac, vector< vector<CpxNumMat> >& c, CpxNumMat& x)
{
  assert(nbscales==c.size() && nbangles_coarse==c[1].size());
  //int F1 = N1/2;  int F2 = N2/2;
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
		fdct_usfft_adjsepangle(XL1, XL2, nbangles[sc], msc, Xhghs[sc]);
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
		fdct_usfft_adjsepangle(XL1, XL2, nbangles[sc], msc, Xhghs[sc]);
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

//------------------------------------
int fdct_usfft_adjfftL(int N1, int N2, CpxOffMat& O, CpxNumMat& x)
{
  int F1 = N1/2;  int F2 = N2/2;
  //ifftshift
  CpxNumMat T(N1, N2);  fdct_usfft_ifftshift(O, T);
  //ifft 
  fftwnd_plan p = fftw2d_create_plan(N2, N1, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fftwnd_one(p, (fftw_complex*)T.data(), NULL);
  fftwnd_destroy_plan(p);
  double sqrtprod = sqrt(double(N1*N2));  for(int i=0; i<N1; i++)	 for(int j=0; j<N2; j++)		T(i,j) /= sqrtprod;
  //move zero freq to front
  x = T; //x.resize(N1, N2, -F1, -F2);  //fdct_usfft_fftshift(N1, N2, F1, F2, T, x);
  return 0;
}

//------------------------------------
int fdct_usfft_adjsepscale(int N1, int N2, int nbscales, int ac, vector<CpxOffMat>& Xhghs, CpxOffMat& O)
{
  int F1 = N1/2;  int F2 = N2/2;
  CpxOffMat X;
  if(ac==1) {
	 double XL1 = 4.0*N1/3.0;  double XL2 = 4.0*N2/3.0; //range
	 int XS1, XS2;  int XF1, XF2;  double XR1, XR2;	 fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
	 X.resize(XS1, XS2);
  } else {
	 X.resize(N1, N2);
  }
  double XL1 = 4.0*N1/3.0;  double XL2 = 4.0*N2/3.0;
  int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
  for(int sc=nbscales-1; sc>0; sc--) {
	 double XL1n = XL1/2;	 double XL2n = XL2/2;
	 int XS1n, XS2n;	 int XF1n, XF2n;	 double XR1n, XR2n;
	 fdct_usfft_rangecompute(XL1n, XL2n, XS1n, XS2n, XF1n, XF2n, XR1n, XR2n);
	 
	 DblOffMat lowpass(XS1n, XS2n);
	 fdct_usfft_lowpasscompute(XL1n, XL2n, lowpass);
	 DblOffMat hghpass(XS1n, XS2n);
	 for(int i=-XF1n; i<-XF1n+XS1n; i++)
		for(int j=-XF2n; j<-XF2n+XS2n; j++)
		  hghpass(i,j) = sqrt(1-lowpass(i,j)*lowpass(i,j));
	 for(int i=-XF1n; i<-XF1n+XS1n; i++)
		for(int j=-XF2n; j<-XF2n+XS2n; j++)
		  Xhghs[sc](i,j) *= hghpass(i,j);
	 for(int i=-XF1n; i<-XF1n+XS1n; i++)
		for(int j=-XF2n; j<-XF2n+XS2n; j++)
		  Xhghs[sc-1](i,j) *= lowpass(i,j);
	 CpxOffMat& G = Xhghs[sc];
	 for(int i=G.s(); i<G.s()+G.m(); i++)
		for(int j=G.t(); j<G.t()+G.n(); j++)
		  X(i,j) += G(i,j);
	 XL1 = XL1/2;	 XL2 = XL2/2;
	 fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
  }
  for(int i=-XF1; i<-XF1+XS1; i++)
	 for(int j=-XF2; j<-XF2+XS2; j++)
		X(i,j) += Xhghs[0](i,j);
  
  //fold if necessary
  O.resize(N1, N2);
  if(ac==1) {
	 double XL1 = 4.0*N1/3.0;  double XL2 = 4.0*N2/3.0;
	 int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
	 //times pou;
	 DblOffMat lowpass(XS1,XS2);
	 fdct_usfft_lowpasscompute(XL1, XL2, lowpass);
	 for(int i=-XF1; i<-XF1+XS1; i++)
		for(int j=-XF2; j<-XF2+XS2; j++)
		  X(i,j) *= lowpass(i,j);
	 IntOffVec t1(XS1);
	 for(int i=-XF1; i<-XF1+XS1; i++)		if(     i<-N1/2) t1(i) = i+int(N1);		else if(i>(N1-1)/2) t1(i) = i-int(N1);		else t1(i) = i;
	 IntOffVec t2(XS2);
	 for(int i=-XF2; i<-XF2+XS2; i++)		if(     i<-N2/2) t2(i) = i+int(N2);		else if(i>(N2-1)/2) t2(i) = i-int(N2);		else t2(i) = i;
	 for(int i=-XF1; i<-XF1+XS1; i++)
		for(int j=-XF2; j<-XF2+XS2; j++)
		  O(t1(i), t2(j)) += X(i,j);
  } else {
	 O = X;
  }
  
  return 0;
}

//------------------------------------
int fdct_usfft_adjsepangle(double XL1, double XL2, int nbangles, vector<CpxOffMat>& msc, CpxOffMat& Xhgh)
{
  int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);
  map<int, fftw_plan> f1map; //forward 1d map, IN_PLACE
  map<int, fftw_plan> b1map; //backward 1d map, IN_PLACE
  
  //allocate space
  Xhgh.resize(XS1, XS2);  clear(Xhgh);
  assert(msc.size()==nbangles);
  int nd = nbangles / 4;
  int nbquadrants = 4;
  int wcnt = 0;
  
  CpxOffMat Xhghb(Xhgh);
  double XL1b = XL1;  double XL2b = XL2;
  
  int qvec[] = {2,1,0,3};
  for(int qi=0; qi<nbquadrants; qi++) {
	 int q = qvec[qi];
	 //rotate
	 fdct_usfft_rotate_forward(q, XL1b, XL2b, XL1, XL2);	 XL1 = abs(XL1);	 XL2 = abs(XL2);
	 fdct_usfft_rotate_forward(q, Xhghb, Xhgh);

	 int XS1, XS2;  int XF1, XF2;  double XR1, XR2;  fdct_usfft_rangecompute(XL1, XL2, XS1, XS2, XF1, XF2, XR1, XR2);	 //double XW1 = XL1/nd;
	 double XW2 = XL2/nd;
	 double xs = XR1/4; 		double xe = XR1;	 int xf = int(ceil(xs));
	 //xn, yn
	 int xn = int(ceil(xe-xs));	 int yn = 2*int(ceil(XW2))+1;
	 if(xn%2==0) xn++;	 if(yn%2==0) yn++;
	 int xh = xn/2;	 int yh = yn/2;
	 vector<CpxOffMat> tsc(nd);
	 //rotate data forward into msc
	 for(int w=nd-1; w>=0; w--) {
		fdct_usfft_rotate_forward(q, msc[wcnt], tsc[w]);
		wcnt++;
	 }
	 //the adjoint of sample using USFFT for each line
	 for(int xcur=xf; xcur<xe; xcur++) { //for each line
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
		//adjoint of interpolation
		int tmpx = xcur%xn;		if(tmpx<-xh) tmpx+=xn;				  if(tmpx>=-xh+xn) tmpx-=xn;
		CpxNumMat res(yn, nd);
		for(int w=0; w<nd; w++)
		  for(int yid=0; yid<yn; yid++) {
			 res(yid,w) = tsc[w](tmpx,yid-yh);
		  }
		CpxOffVec val(XS2);
		fdct_usfft_adj1dinterp(ysft, wsft, res, val, f1map, b1map);
		for(int i=-XF2; i<-XF2+XS2; i++) {
		  Xhgh(xcur,i) += val(i);
		}
	 }
	 fdct_usfft_rotate_backward(q, Xhgh, Xhghb);
  } //quadrant
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

//------------------------------------
int fdct_usfft_adjifftS(vector<CpxNumMat>& csc, vector<CpxOffMat>& msc)
{
  typedef pair<int,int> intpair;
  map<intpair, fftwnd_plan> planmap;
  //do work
  msc.resize(csc.size());
  for(int w=0; w<csc.size(); w++) {
	 int xn = csc[w].m();	 int yn = csc[w].n();
	 int xh = xn/2;	 int yh = yn/2;

	 CpxNumMat tpdata = csc[w];	 //shift	 //CpxNumMat tpdata(xn,yn);	 fdct_usfft_ifftshift(xn,yn,xh,yh,csc[w],tpdata);
	 //fft
	 map<intpair,fftwnd_plan>::iterator mit=planmap.find( intpair(xn,yn) );
	 fftwnd_plan p = NULL;
	 if(mit!=planmap.end()) {
		p = (*mit).second;
	 } else {
		p = fftw2d_create_plan(yn, xn, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
		planmap[ intpair(xn,yn) ] = p;
	 }
	 fftwnd_one(p, (fftw_complex*)tpdata.data(), NULL);
	 double sqrtprod = sqrt(double(xn*yn));
	 for(int i=0; i<xn; i++)		for(int j=0; j<yn; j++)		  tpdata(i,j) /= sqrtprod;
	 //shift
	 msc[w].resize(xn,yn);	 fdct_usfft_fftshift(tpdata,msc[w]);
  }
  for(map<intpair,fftwnd_plan>::iterator mit=planmap.begin(); mit!=planmap.end(); mit++) {
	 fftwnd_plan p = (*mit).second;
	 fftwnd_destroy_plan(p);
  }
  return 0;
}

//------------------------------------
int fdct_usfft_adjwavelet(vector<CpxNumMat>& csc, CpxOffMat& Xhgh)
{
  assert(csc.size()==1);
  CpxNumMat& C = csc[0];
  int N1 = C.m();  int N2 = C.n();
  
  CpxNumMat T(C);  //ifftshift  //CpxNumMat T(N1, N2);  fdct_usfft_ifftshift(N1, N2, F1, F2, csc[0], T);
  //fft
  fftwnd_plan p = fftw2d_create_plan(N2, N1, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fftwnd_one(p, (fftw_complex*)T.data(), NULL);
  fftwnd_destroy_plan(p);
  double sqrtprod = sqrt(double(N1*N2));  for(int i=0; i<N1; i++)	 for(int j=0; j<N2; j++)		T(i,j) /= sqrtprod;
  Xhgh.resize(N1, N2);
  fdct_usfft_fftshift(T, Xhgh);
  return 0;
}

//------------------------------------
int fdct_usfft_adj1dinterp(DblNumMat& off, DblNumMat& wgt, CpxNumMat& res, CpxOffVec& val,
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
	 map<int, fftw_plan>::iterator bit = b1map.find(N);
	 fftw_plan bp = NULL;
	 if(bit!=b1map.end()) {	 bp = (*bit).second;
	 } else {	 bp = fftw_create_plan(N, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);	 b1map[N] = bp;
	 }
	 CpxOffVec feq(N);
	 CpxOffVec vvv(N);
	 CpxOffVec fff(N);
	 CpxNumVec tmp(N);
	 //1. for each wedge, from res, fft, multiply, add back
	 int nbwedges = off.n();
	 CpxOffVec sum(N); //sum, equal to zero
	 for(int w=0; w<nbwedges; w++) {
		for(int k=0; k<off.m(); k++) {
		  int rk = (k>F)? k-N : k;
		  vvv(rk) = res(k,w) * wgt(k,w);
		} //new value
		fdct_usfft_ifftshift(vvv, tmp);
		fftw_one(fp, (fftw_complex*)tmp.data(), NULL);
		double scale = sqrt(double(N));	 for(int k=0; k<N; k++)		  tmp(k) /= scale;	  //scaling
		fdct_usfft_fftshift(tmp, fff); ///new freq
		for(int k=-F; k<-F+N; k++) {
		  double phase = -2*M_PI*k/double(N) * off(0,w);
		  feq(k) = fff(k) * polar(1.0, phase);
		} //freq
		for(int k=-F; k<-F+N; k++) {
		  sum(k) += feq(k);
		}
	 }
	 //2. fft_mid val
	 fdct_usfft_ifftshift(sum, tmp);
	 fftw_one(bp, (fftw_complex*)tmp.data(), NULL);
	 double scale = sqrt(double(N));  for(int k=0; k<N; k++) {	 tmp(k) /= scale;	 } //scaling
	 fdct_usfft_fftshift(tmp, val);
  } else { //USFFT INTERPOLATION
  	 //--------------------------------------------
	 int L = 4;
	 int D = 16;
	 int N = val.m();
	 int F = -val.s();
	 fftw_plan fp = NULL; //length D*N
	 map<int, fftw_plan>::iterator fit=f1map.find(D*N);
	 if(fit!=f1map.end()) {		fp = (*fit).second;
	 } else {	 fp = fftw_create_plan(D*N, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);	 f1map[D*N] = fp;
	 }
	 fftw_plan bp = NULL; //length N
	 map<int, fftw_plan>::iterator bit = b1map.find(N);
	 if(bit!=b1map.end()) {	 bp = (*bit).second;
	 } else {	 bp = fftw_create_plan(N, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);	 b1map[N] = bp;
	 }
	 CpxNumVec tmp(N);
	 CpxNumVec exttmp(D*N);
	 CpxOffVec feq(N);	 clear(feq);
	 CpxOffVec pok(N);	 setvalue(pok, cpx(1,0)); //power of k, unit vector
	 vector<CpxOffVec> extdat(L);	 //CpxOffVec extdat[4];
	 for(int l=0; l<L; l++)	{
		extdat[l].resize(D*N);
	 }
	 //4. back taylor
	 double step = 2*M_PI/(D*N);
	 for(int j=0; j<off.n(); j++)
		for(int i=0; i<off.m(); i++) {
		  double cof = off(i,j) / N * 2*M_PI;
		  if(cof<-M_PI)			 cof+= 2*M_PI;		  else if(cof>=M_PI)			 cof-= 2*M_PI; //collapse into [-pi,pi)
		  int ind = int(floor(cof/step));		  ind = min(max(ind, -D*N/2), D*N/2-1); 		  //assert(ind>=-D*N/2 && ind<D*N/2);
		  double dta = cof-ind*step;
		  double pow = 1.0;
		  cpx sss( res(i,j) * wgt(i,j) );
		  for(int l=0; l<L; l++) {
			 extdat[l](ind) += sss * pow;
			 pow = pow * dta/(l+1);
		  }
		}
	 for(int l=0; l<L; l++) {
		fdct_usfft_ifftshift(extdat[l], exttmp);
		fftw_one(fp, (fftw_complex*)exttmp.data(), NULL);
		fdct_usfft_fftshift(exttmp, extdat[l]);
	 }
	 for(int l=0; l<L; l++) {		//1. sum		//2. update pok
		for(int k=-F; k<-F+N; k++) {		  feq(k) += extdat[l](k) * pok(k);		}
		for(int k=-F; k<-F+N; k++) {		  pok(k) *= cpx(0,-k);		}
	 }
	 fdct_usfft_ifftshift(feq, tmp);
	 fftw_one(bp, (fftw_complex*)tmp.data(), NULL);
	 double scale = double(N);	 for(int k=0; k<N; k++) {	 tmp(k) /= scale;	 } //scaling
	 fdct_usfft_fftshift(tmp, val);
  }
  return 0;
}

FDCT_USFFT_NS_END_NAMESPACE
