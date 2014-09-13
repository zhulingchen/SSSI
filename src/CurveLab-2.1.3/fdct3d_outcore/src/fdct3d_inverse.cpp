/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#include "fdct3d.hpp"
#include "fdct3dinline.hpp"

int fdct3d_inverse_angles(double L1,double L2,double L3, int s,int nd, CpxCrvletOcr& C, CpxOffTns& O)
{
  //vector<CpxNumTns>& csc = C[s];
  
  int nf = 6;
  int wcnt = 0;
  int S1, S2, S3;	 int F1, F2, F3;	 double R1, R2, R3;	 fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);
  DblOffVec big1(S1);  fdct3d_lowpass(L1, big1);
  DblOffVec big2(S2);  fdct3d_lowpass(L2, big2);
  DblOffVec big3(S3);  fdct3d_lowpass(L3, big3);
  
  double Lh1 = L1/2;  double Lh2 = L2/2;  double Lh3 = L3/2;
  int Sh1, Sh2, Sh3;	 int Fh1, Fh2, Fh3;	 double Rh1, Rh2, Rh3;	 fdct3d_rangecompute(Lh1, Lh2, Lh3, Sh1, Sh2, Sh3, Fh1, Fh2, Fh3, Rh1, Rh2, Rh3);
  DblOffVec sma1(S1);  fdct3d_lowpass(Lh1, sma1);
  DblOffVec sma2(S2);  fdct3d_lowpass(Lh2, sma2);
  DblOffVec sma3(S3);  fdct3d_lowpass(Lh3, sma3);
  
  double W1 = L1/nd;  double W2 = L2/nd;  double W3 = L3/nd;
  
  typedef pair<int,int> intpair;  typedef pair<int, intpair> inttriple;
  map<inttriple, fftwnd_plan> planmap;

  //face 0: x,y,z
  for(int h=0; h<nd; h++) { //(y first z second)
	 for(int g=0; g<nd; g++) {

		  double xs = R1/4-(W1/2)/4;		double xe = R1;
		  double ys = -R2 + (2*g-1)*W2/2;		double ye = -R2 + (2*g+3)*W2/2;
		  double zs = -R3 + (2*h-1)*W3/2;		double ze = -R3 + (2*h+3)*W3/2;
		  int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		  double thts, thtm, thte; //y to x
		  if(g==0) {
			 thts = atan2(-1.0, 1.0-1.0/nd);			 thtm = atan2(-1.0+1.0/nd, 1.0);			 thte = atan2(-1.0+3.0/nd, 1.0);
		  } else if(g==nd-1) {
			 thts = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);			 thtm = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);			 thte = atan2(1.0, 1.0-1.0/nd);
		  } else {
			 thts = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);			 thtm = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);			 thte = atan2(-1.0+(2.0*g+3.0)/nd, 1.0);
		  }
		  double phis, phim, phie; //z to x
		  if(h==0) {
			 phis = atan2(-1.0, 1.0-1.0/nd);			 phim = atan2(-1.0+1.0/nd, 1.0);			 phie = atan2(-1.0+3.0/nd, 1.0);
		  } else if(h==nd-1) {
			 phis = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);			 phim = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);			 phie = atan2(1.0, 1.0-1.0/nd);
		  } else {
			 phis = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);			 phim = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);			 phie = atan2(-1.0+(2.0*h+3.0)/nd, 1.0);
		  }
		  int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2; //half
		  double R21 = R2/R1;		  double R31 = R3/R1;
		  
		  CpxNumTns tpdata(xn,yn,zn);
		  //CpxNumTns& Cblk = C.block(s,wcnt);
		  //tpdata = Cblk;
		  ///tpdata = csc[wcnt];
		  tpdata = C.block(s,wcnt);
		  //fft
		  fftwnd_plan p = NULL;
		  map<inttriple, fftwnd_plan>::iterator mit = planmap.find( inttriple(xn, intpair(yn,zn)) );
		  if(mit!=planmap.end()) {			 p = (*mit).second;
		  } else {
			 p = fftw3d_create_plan(zn, yn, xn, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
			 planmap[ inttriple(xn, intpair(yn,zn)) ] = p;
		  }
		  fftwnd_one(p, (fftw_complex*)tpdata.data(), NULL);		  //cerr<<"wedge s"<<endl;
		  double sqrtprod = sqrt(double(xn*yn*zn));
		  for(int i=0; i<xn; i++)			 for(int j=0; j<yn; j++)				for(int k=0; k<zn; k++)				  tpdata(i,j,k) /= sqrtprod;
		  CpxOffTns wpdata(xn,yn,zn);
		  fdct3d_fftshift(xn,yn,zn,tpdata,wpdata);
		  
		  for(int xcur=(int)ceil(xs); xcur<xe; xcur++) {
			 int yfm = (int)ceil( max(-R2, R21*xcur*tan(thts)) );
			 int yto = (int)floor( min(R2, R21*xcur*tan(thte)) );
			 int zfm = (int)ceil( max(-R3, R31*xcur*tan(phis)) );
			 int zto = (int)floor( min(R3, R31*xcur*tan(phie)) );
			 for(int ycur=yfm; ycur<=yto; ycur++)
				for(int zcur=zfm; zcur<=zto; zcur++) {
				  int tmpx = xcur%xn;				  if(tmpx<-xh) tmpx+=xn;				  if(tmpx>=-xh+xn) tmpx-=xn;
				  int tmpy = ycur%yn;				  if(tmpy<-yh) tmpy+=yn;				  if(tmpy>=-yh+yn) tmpy-=yn;
				  int tmpz = zcur%zn;				  if(tmpz<-zh) tmpz+=zn;				  if(tmpz>=-zh+zn) tmpz-=zn;
				  
				  double thtcur = atan2(ycur/R2, xcur/R1);
				  double phicur = atan2(zcur/R3, xcur/R1);
				  double glbpou;					 fdct3d_globalpou(thtcur, phicur, M_PI/4-atan2(1.0-1.0/nd, 1.0), glbpou);
				  double wtht;
				  if(thtcur<thtm) {
					 if(g==0)						wtht = 1;
					 else {						double l,r;						fdct3d_window( (thtcur-thts)/(thtm-thts), l, r);						wtht = l;					 }
				  } else {
					 if(g==nd-1)						wtht = 1;
					 else {						double l,r;						  fdct3d_window( (thtcur-thtm)/(thte-thtm), l, r);						wtht = r;					 }
				  }
				  double wphi;
				  if(phicur<phim) {
					 if(h==0)						wphi = 1;
					 else {						double l,r;						  fdct3d_window( (phicur-phis)/(phim-phis), l, r);						wphi = l;					 }
				  } else {
					 if(h==nd-1)						wphi = 1;
					 else {						double l,r;						  fdct3d_window( (phicur-phim)/(phie-phim), l, r);						wphi = r;					 }
				  }
				  double pou = glbpou * wtht * wphi;
				  wpdata(tmpx, tmpy, tmpz) *= pou;
				  
				  double ss = sma1(xcur)*sma2(ycur)*sma3(zcur);				  double bb = big1(xcur)*big2(ycur)*big3(zcur);
				  //int bi,bj,bk;			 int oi,oj,ok;			 fdct3d_position_aux(N1,N2,N3,b, xcur,ycur,zcur, bi,bj,bk,oi,oj,ok);
				  //CpxNumTns& Wblk = W.block(bi,bj,bk);
				  //Wblk(oi,oj,ok) += wpdata(tmpx,tmpy,tmpz)  * bb * sqrt(1.0-ss*ss);
				  O(xcur,ycur,zcur) += wpdata(tmpx,tmpy,tmpz)  * bb * sqrt(1.0-ss*ss);
				}
		  } //xcur

		wcnt++;
	 }
  } //end of face
  //face 1. y z x
  for(int f=0; f<nd; f++) {
	 for(int h=0; h<nd; h++) {

		  double ys = R2/4-(W2/2)/4;		  double ye = R2;
		  double zs = -R3 + (2*h-1)*W3/2;		  double ze = -R3 + (2*h+3)*W3/2;
		  double xs = -R1 + (2*f-1)*W1/2;		  double xe = -R1 + (2*f+3)*W1/2;
		  int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		  double thts, thtm, thte; //z to y
		  if(h==0) {
			 thts = atan2(-1.0, 1.0-1.0/nd);			 thtm = atan2(-1.0+1.0/nd, 1.0);			 thte = atan2(-1.0+3.0/nd, 1.0);
		  } else if(h==nd-1) {
			 thts = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);			 thtm = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);			 thte = atan2(1.0, 1.0-1.0/nd);
		  } else {
			 thts = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);			 thtm = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);			 thte = atan2(-1.0+(2.0*h+3.0)/nd, 1.0);
		  }
		  double phis, phim, phie; //z to x
		  if(f==0) {
			 phis = atan2(-1.0, 1.0-1.0/nd);			 phim = atan2(-1.0+1.0/nd, 1.0);			 phie = atan2(-1.0+3.0/nd, 1.0);
		  } else if(f==nd-1) {
			 phis = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);			 phim = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);			 phie = atan2(1.0, 1.0-1.0/nd);
		  } else {
			 phis = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);			 phim = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);			 phie = atan2(-1.0+(2.0*f+3.0)/nd, 1.0);
		  }
		  int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2;
		  double R32 = R3/R2;		  double R12 = R1/R2;
		  
		  CpxNumTns tpdata(xn,yn,zn);
		  //CpxNumTns& Cblk = C.block(s,wcnt);
		  //tpdata = Cblk;
		  ///tpdata = csc[wcnt];
		  tpdata = C.block(s,wcnt);
		  //fft
		  fftwnd_plan p = NULL;
		  map<inttriple, fftwnd_plan>::iterator mit = planmap.find( inttriple(xn, intpair(yn,zn)) );
		  if(mit!=planmap.end()) {			 p = (*mit).second;
		  } else {
			 p = fftw3d_create_plan(zn, yn, xn, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
			 planmap[ inttriple(xn, intpair(yn,zn)) ] = p;
		  }
		  fftwnd_one(p, (fftw_complex*)tpdata.data(), NULL);		  //cerr<<"wedge s"<<endl;
		  double sqrtprod = sqrt(double(xn*yn*zn));
		  for(int i=0; i<xn; i++)			 for(int j=0; j<yn; j++)				for(int k=0; k<zn; k++)				  tpdata(i,j,k) /= sqrtprod;
		  CpxOffTns wpdata(xn,yn,zn);
		  fdct3d_fftshift(xn,yn,zn,tpdata,wpdata);

		  for(int ycur=(int)ceil(ys); ycur<ye; ycur++) {
			 int zfm = (int)ceil( max(-R3, R32*ycur*tan(thts)) );
			 int zto = (int)floor( min(R3, R32*ycur*tan(thte)) );
			 int xfm = (int)ceil( max(-R1, R12*ycur*tan(phis)) );
			 int xto = (int)floor( min(R1, R12*ycur*tan(phie)) );
			 for(int zcur=zfm; zcur<=zto; zcur++)
				for(int xcur=xfm; xcur<=xto; xcur++) {
				  int tmpx = xcur%xn;				  if(tmpx<-xh) tmpx+=xn;				  if(tmpx>=-xh+xn) tmpx-=xn;
				  int tmpy = ycur%yn;				  if(tmpy<-yh) tmpy+=yn;				  if(tmpy>=-yh+yn) tmpy-=yn;
				  int tmpz = zcur%zn;				  if(tmpz<-zh) tmpz+=zn;				  if(tmpz>=-zh+zn) tmpz-=zn;
				  
				  double thtcur = atan2(zcur/R3, ycur/R2);
				  double phicur = atan2(xcur/R1, ycur/R2);
				  double glbpou;					 fdct3d_globalpou(thtcur, phicur, M_PI/4-atan2(1.0-1.0/nd, 1.0), glbpou); //CHECK
				  double wtht;
				  if(thtcur<thtm) {
					 if(h==0)						wtht = 1;
					 else {						double l,r;						fdct3d_window( (thtcur-thts)/(thtm-thts), l, r);						wtht = l;					 }
				  } else {
					 if(h==nd-1)						wtht = 1;
					 else {						double l,r;						  fdct3d_window( (thtcur-thtm)/(thte-thtm), l, r);						wtht = r;					 }
				  }
				  double wphi;
				  if(phicur<phim) {
					 if(f==0)						wphi = 1;
					 else {						double l,r;						  fdct3d_window( (phicur-phis)/(phim-phis), l, r);						wphi = l;					 }
				  } else {
					 if(f==nd-1)						wphi = 1;
					 else {						double l,r;						  fdct3d_window( (phicur-phim)/(phie-phim), l, r);						wphi = r;					 }
				  }
				  double pou = glbpou * wtht * wphi;
				  wpdata(tmpx, tmpy, tmpz) *= pou;
				  
				  double ss = sma1(xcur)*sma2(ycur)*sma3(zcur);				  double bb = big1(xcur)*big2(ycur)*big3(zcur);
				  //int bi,bj,bk;			 int oi,oj,ok;			 fdct3d_position_aux(N1,N2,N3,b, xcur,ycur,zcur, bi,bj,bk,oi,oj,ok);
				  //CpxNumTns& Wblk = W.block(bi,bj,bk);
				  //Wblk(oi,oj,ok) += wpdata(tmpx,tmpy,tmpz)  * bb * sqrt(1.0-ss*ss);
				  O(xcur,ycur,zcur) += wpdata(tmpx,tmpy,tmpz)  * bb * sqrt(1.0-ss*ss);
				}
		  } //ycur

		wcnt++;
	 }
  }//end of face
  //face 2. z x y
  for(int g=0; g<nd; g++) {
	 for(int f=0; f<nd; f++) {
		
		  double zs = R3/4-(W3/2)/4;		double ze = R3;
		  double xs = -R1 + (2*f-1)*W1/2;		double xe = -R1 + (2*f+3)*W1/2;
		  double ys = -R2 + (2*g-1)*W2/2;		double ye = -R2 + (2*g+3)*W2/2;
		  int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		  double thts, thtm, thte; //y to x
		  if(f==0) {
			 thts = atan2(-1.0, 1.0-1.0/nd);			 thtm = atan2(-1.0+1.0/nd, 1.0);			 thte = atan2(-1.0+3.0/nd, 1.0);
		  } else if(f==nd-1) {
			 thts = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);			 thtm = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);			 thte = atan2(1.0, 1.0-1.0/nd);
		  } else {
			 thts = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);			 thtm = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);			 thte = atan2(-1.0+(2.0*f+3.0)/nd, 1.0);
		  }
		  double phis, phim, phie; //z to x
		  if(g==0) {
			 phis = atan2(-1.0, 1.0-1.0/nd);			 phim = atan2(-1.0+1.0/nd, 1.0);			 phie = atan2(-1.0+3.0/nd, 1.0);
		  } else if(g==nd-1) {
			 phis = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);			 phim = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);			 phie = atan2(1.0, 1.0-1.0/nd);
		  } else {
			 phis = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);			 phim = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);			 phie = atan2(-1.0+(2.0*g+3.0)/nd, 1.0);
		  }
		  int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2;
		  double R13 = R1/R3;		  double R23 = R2/R3;//double R13 = double(F1)/double(F3);		  double R23 = double(F2)/double(F3);
		  
		  CpxNumTns tpdata(xn,yn,zn);
		  //CpxNumTns& Cblk = C.block(s,wcnt);
		  //tpdata = Cblk;
		  ///tpdata = csc[wcnt];
		  tpdata = C.block(s,wcnt);
		  //fft
		  fftwnd_plan p = NULL;
		  map<inttriple, fftwnd_plan>::iterator mit = planmap.find( inttriple(xn, intpair(yn,zn)) );
		  if(mit!=planmap.end()) {			 p = (*mit).second;
		  } else {
			 p = fftw3d_create_plan(zn, yn, xn, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
			 planmap[ inttriple(xn, intpair(yn,zn)) ] = p;
		  }
		  fftwnd_one(p, (fftw_complex*)tpdata.data(), NULL);		  //cerr<<"wedge s"<<endl;
		  double sqrtprod = sqrt(double(xn*yn*zn));
		  for(int i=0; i<xn; i++)			 for(int j=0; j<yn; j++)				for(int k=0; k<zn; k++)				  tpdata(i,j,k) /= sqrtprod;
		  CpxOffTns wpdata(xn,yn,zn);
		  fdct3d_fftshift(xn,yn,zn,tpdata,wpdata);
		  
		  for(int zcur=(int)ceil(zs); zcur<ze; zcur++) {
			 int xfm = (int)ceil( max(-R1, R13*zcur*tan(thts)) );
			 int xto = (int)floor( min(R1, R13*zcur*tan(thte)) );
			 int yfm = (int)ceil( max(-R2, R23*zcur*tan(phis)) );
			 int yto = (int)floor( min(R2, R23*zcur*tan(phie)) );
			 for(int xcur=xfm; xcur<=xto; xcur++)
				for(int ycur=yfm; ycur<=yto; ycur++) {
				  int tmpx = xcur%xn;				  if(tmpx<-xh) tmpx+=xn;				  if(tmpx>=-xh+xn) tmpx-=xn;
				  int tmpy = ycur%yn;				  if(tmpy<-yh) tmpy+=yn;				  if(tmpy>=-yh+yn) tmpy-=yn;
				  int tmpz = zcur%zn;				  if(tmpz<-zh) tmpz+=zn;				  if(tmpz>=-zh+zn) tmpz-=zn;
				  
				  double thtcur = atan2(xcur/R1, zcur/R3);
				  double phicur = atan2(ycur/R2, zcur/R3);
				  double glbpou;					 fdct3d_globalpou(thtcur, phicur, M_PI/4-atan2(1.0-1.0/nd, 1.0), glbpou);
				  double wtht;
				  if(thtcur<thtm) {
					 if(f==0)						wtht = 1;
					 else {						double l,r;						fdct3d_window( (thtcur-thts)/(thtm-thts), l, r);						wtht = l;					 }
				  } else {
					 if(f==nd-1)						wtht = 1;
					 else {						double l,r;						  fdct3d_window( (thtcur-thtm)/(thte-thtm), l, r);						wtht = r;					 }
				  }
				  double wphi;
				  if(phicur<phim) {
					 if(g==0)						wphi = 1;
					 else {						double l,r;						  fdct3d_window( (phicur-phis)/(phim-phis), l, r);						wphi = l;					 }
				  } else {
					 if(g==nd-1)						wphi = 1;
					 else {						double l,r;						  fdct3d_window( (phicur-phim)/(phie-phim), l, r);						wphi = r;					 }
				  }
				  double pou = glbpou * wtht * wphi;
				  wpdata(tmpx, tmpy, tmpz) *= pou;
				  
				  double ss = sma1(xcur)*sma2(ycur)*sma3(zcur);				  double bb = big1(xcur)*big2(ycur)*big3(zcur);
				  //int bi,bj,bk;			 int oi,oj,ok;			 fdct3d_position_aux(N1,N2,N3,b, xcur,ycur,zcur, bi,bj,bk,oi,oj,ok);
				  //CpxNumTns& Wblk = W.block(bi,bj,bk);
				  //Wblk(oi,oj,ok) += wpdata(tmpx,tmpy,tmpz)  * bb * sqrt(1.0-ss*ss);
				  O(xcur,ycur,zcur) += wpdata(tmpx,tmpy,tmpz)  * bb * sqrt(1.0-ss*ss);
				}
		  }//zcur

		wcnt++;
	 }
  }//end of face
  //face 3: -x,-y,-z
  for(int h=nd-1; h>=0; h--) {
	 for(int g=nd-1; g>=0; g--) {

		  double xs = -R1;		  double xe = -R1/4+(W1/2)/4;
		  double ys = -R2 + (2*g-1)*W2/2;		double ye = -R2 + (2*g+3)*W2/2;
		  double zs = -R3 + (2*h-1)*W3/2;		double ze = -R3 + (2*h+3)*W3/2;
		  int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		  double thts, thtm, thte; //y to x
		  if(g==0) {
			 thts = atan2(-1.0, 1.0-1.0/nd);			 thtm = atan2(-1.0+1.0/nd, 1.0);			 thte = atan2(-1.0+3.0/nd, 1.0);
		  } else if(g==nd-1) {
			 thts = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);			 thtm = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);			 thte = atan2(1.0, 1.0-1.0/nd);
		  } else {
			 thts = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);			 thtm = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);			 thte = atan2(-1.0+(2.0*g+3.0)/nd, 1.0);
		  }
		  double phis, phim, phie; //z to x
		  if(h==0) {
			 phis = atan2(-1.0, 1.0-1.0/nd);			 phim = atan2(-1.0+1.0/nd, 1.0);			 phie = atan2(-1.0+3.0/nd, 1.0);
		  } else if(h==nd-1) {
			 phis = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);			 phim = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);			 phie = atan2(1.0, 1.0-1.0/nd);
		  } else {
			 phis = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);			 phim = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);			 phie = atan2(-1.0+(2.0*h+3.0)/nd, 1.0);
		  }
		  int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2;
		  double R21 = R2/R1;		  double R31 = R3/R1;
		  
		  CpxNumTns tpdata(xn,yn,zn);
		  //CpxNumTns& Cblk = C.block(s,wcnt);
		  //tpdata = Cblk;
		  ///tpdata = csc[wcnt];
		  tpdata = C.block(s,wcnt);
		  //fft
		  fftwnd_plan p = NULL;
		  map<inttriple, fftwnd_plan>::iterator mit = planmap.find( inttriple(xn, intpair(yn,zn)) );
		  if(mit!=planmap.end()) {			 p = (*mit).second;
		  } else {
			 p = fftw3d_create_plan(zn, yn, xn, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
			 planmap[ inttriple(xn, intpair(yn,zn)) ] = p;
		  }
		  fftwnd_one(p, (fftw_complex*)tpdata.data(), NULL);		  //cerr<<"wedge s"<<endl;
		  double sqrtprod = sqrt(double(xn*yn*zn));
		  for(int i=0; i<xn; i++)			 for(int j=0; j<yn; j++)				for(int k=0; k<zn; k++)				  tpdata(i,j,k) /= sqrtprod;
		  CpxOffTns wpdata(xn,yn,zn);
		  fdct3d_fftshift(xn,yn,zn,tpdata,wpdata);
		  
		  for(int xcur=(int)ceil(xs); xcur<xe; xcur++) {
			 int yfm = (int)ceil( max(-R2, R21*(-xcur)*tan(thts)) );
			 int yto = (int)floor( min(R2, R21*(-xcur)*tan(thte)) );
			 int zfm = (int)ceil( max(-R3, R31*(-xcur)*tan(phis)) );
			 int zto = (int)floor( min(R3, R31*(-xcur)*tan(phie)) );
			 for(int ycur=yfm; ycur<=yto; ycur++)
				for(int zcur=zfm; zcur<=zto; zcur++) {
				  int tmpx = xcur%xn;				  if(tmpx<-xh) tmpx+=xn;				  if(tmpx>=-xh+xn) tmpx-=xn;
				  int tmpy = ycur%yn;				  if(tmpy<-yh) tmpy+=yn;				  if(tmpy>=-yh+yn) tmpy-=yn;
				  int tmpz = zcur%zn;				  if(tmpz<-zh) tmpz+=zn;				  if(tmpz>=-zh+zn) tmpz-=zn;
				  
				  double thtcur = atan2(ycur/R2, (-xcur)/R1);
				  double phicur = atan2(zcur/R3, (-xcur)/R1);
				  double glbpou;					 fdct3d_globalpou(thtcur, phicur, M_PI/4-atan2(1.0-1.0/nd, 1.0), glbpou);
				  double wtht;
				  if(thtcur<thtm) {
					 if(g==0)						wtht = 1;
					 else {						double l,r;						fdct3d_window( (thtcur-thts)/(thtm-thts), l, r);						wtht = l;					 }
				  } else {
					 if(g==nd-1)						wtht = 1;
					 else {						double l,r;						  fdct3d_window( (thtcur-thtm)/(thte-thtm), l, r);						wtht = r;					 }
				  }
				  double wphi;
				  if(phicur<phim) {
					 if(h==0)						wphi = 1;
					 else {						double l,r;						  fdct3d_window( (phicur-phis)/(phim-phis), l, r);						wphi = l;					 }
				  } else {
					 if(h==nd-1)						wphi = 1;
					 else {						double l,r;						  fdct3d_window( (phicur-phim)/(phie-phim), l, r);						wphi = r;					 }
				  }
				  double pou = glbpou * wtht * wphi;
				  wpdata(tmpx, tmpy, tmpz) *= pou;
				  
				  double ss = sma1(xcur)*sma2(ycur)*sma3(zcur);				  double bb = big1(xcur)*big2(ycur)*big3(zcur);
				  //int bi,bj,bk;			 int oi,oj,ok;			 fdct3d_position_aux(N1,N2,N3,b, xcur,ycur,zcur, bi,bj,bk,oi,oj,ok);
				  //CpxNumTns& Wblk = W.block(bi,bj,bk);
				  //Wblk(oi,oj,ok) += wpdata(tmpx,tmpy,tmpz)  * bb * sqrt(1.0-ss*ss);
				  O(xcur,ycur,zcur) += wpdata(tmpx,tmpy,tmpz)  * bb * sqrt(1.0-ss*ss);
				}
		  } //xcur

		wcnt++;
	 }
  } //end of face
  //face 4: -y,-z,-x
  for(int f=nd-1; f>=0; f--) {
	 for(int h=nd-1; h>=0; h--) {

		  double ys = -R2;		  double ye = -R2/4+(W2/2)/4;
		  double zs = -R3 + (2*h-1)*W3/2;		  double ze = -R3 + (2*h+3)*W3/2;
		  double xs = -R1 + (2*f-1)*W1/2;		  double xe = -R1 + (2*f+3)*W1/2;
		  int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		  double thts, thtm, thte; //z to y
		  if(h==0) {
			 thts = atan2(-1.0, 1.0-1.0/nd);			 thtm = atan2(-1.0+1.0/nd, 1.0);			 thte = atan2(-1.0+3.0/nd, 1.0);
		  } else if(h==nd-1) {
			 thts = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);			 thtm = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);			 thte = atan2(1.0, 1.0-1.0/nd);
		  } else {
			 thts = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);			 thtm = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);			 thte = atan2(-1.0+(2.0*h+3.0)/nd, 1.0);
		  }
		  double phis, phim, phie; //z to x
		  if(f==0) {
			 phis = atan2(-1.0, 1.0-1.0/nd);			 phim = atan2(-1.0+1.0/nd, 1.0);			 phie = atan2(-1.0+3.0/nd, 1.0);
		  } else if(f==nd-1) {
			 phis = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);			 phim = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);			 phie = atan2(1.0, 1.0-1.0/nd);
		  } else {
			 phis = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);			 phim = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);			 phie = atan2(-1.0+(2.0*f+3.0)/nd, 1.0);
		  }
		  int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2;
		  double R32 = R3/R2;		  double R12 = R1/R2;	  //double R32 = double(F3)/double(F2);		  double R12 = double(F1)/double(F2);
		  
		  CpxNumTns tpdata(xn,yn,zn);
		  //CpxNumTns& Cblk = C.block(s,wcnt);
		  //tpdata = Cblk;
		  ///tpdata = csc[wcnt];
		  tpdata = C.block(s,wcnt);
		  //fft
		  fftwnd_plan p = NULL;
		  map<inttriple, fftwnd_plan>::iterator mit = planmap.find( inttriple(xn, intpair(yn,zn)) );
		  if(mit!=planmap.end()) {			 p = (*mit).second;
		  } else {
			 p = fftw3d_create_plan(zn, yn, xn, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
			 planmap[ inttriple(xn, intpair(yn,zn)) ] = p;
		  }
		  fftwnd_one(p, (fftw_complex*)tpdata.data(), NULL);		  //cerr<<"wedge s"<<endl;
		  double sqrtprod = sqrt(double(xn*yn*zn));
		  for(int i=0; i<xn; i++)			 for(int j=0; j<yn; j++)				for(int k=0; k<zn; k++)				  tpdata(i,j,k) /= sqrtprod;
		  CpxOffTns wpdata(xn,yn,zn);
		  fdct3d_fftshift(xn,yn,zn,tpdata,wpdata);
		  
		  for(int ycur=(int)ceil(ys); ycur<ye; ycur++) {
			 int zfm = (int)ceil( max(-R3, R32*(-ycur)*tan(thts)) );
			 int zto = (int)floor( min(R3, R32*(-ycur)*tan(thte)) );
			 int xfm = (int)ceil( max(-R1, R12*(-ycur)*tan(phis)) );
			 int xto = (int)floor( min(R1, R12*(-ycur)*tan(phie)) );
			 for(int zcur=zfm; zcur<=zto; zcur++)
				for(int xcur=xfm; xcur<=xto; xcur++) {
				  int tmpx = xcur%xn;				  if(tmpx<-xh) tmpx+=xn;				  if(tmpx>=-xh+xn) tmpx-=xn;
				  int tmpy = ycur%yn;				  if(tmpy<-yh) tmpy+=yn;				  if(tmpy>=-yh+yn) tmpy-=yn;
				  int tmpz = zcur%zn;				  if(tmpz<-zh) tmpz+=zn;				  if(tmpz>=-zh+zn) tmpz-=zn;
				  
				  double thtcur = atan2(zcur/R3, (-ycur)/R2);
				  double phicur = atan2(xcur/R1, (-ycur)/R2);
				  double glbpou;					 fdct3d_globalpou(thtcur, phicur, M_PI/4-atan2(1.0-1.0/nd, 1.0), glbpou); //CHECK
				  double wtht;
				  if(thtcur<thtm) {
					 if(h==0)						wtht = 1;
					 else {						double l,r;						fdct3d_window( (thtcur-thts)/(thtm-thts), l, r);						wtht = l;					 }
				  } else {
					 if(h==nd-1)						wtht = 1;
					 else {						double l,r;						  fdct3d_window( (thtcur-thtm)/(thte-thtm), l, r);						wtht = r;					 }
				  }
				  double wphi;
				  if(phicur<phim) {
					 if(f==0)						wphi = 1;
					 else {						double l,r;						  fdct3d_window( (phicur-phis)/(phim-phis), l, r);						wphi = l;					 }
				  } else {
					 if(f==nd-1)						wphi = 1;
					 else {						double l,r;						  fdct3d_window( (phicur-phim)/(phie-phim), l, r);						wphi = r;					 }
				  }
				  double pou = glbpou * wtht * wphi;
				  wpdata(tmpx, tmpy, tmpz) *= pou;
				  
				  double ss = sma1(xcur)*sma2(ycur)*sma3(zcur);				  double bb = big1(xcur)*big2(ycur)*big3(zcur);
				  //int bi,bj,bk;			 int oi,oj,ok;			 fdct3d_position_aux(N1,N2,N3,b, xcur,ycur,zcur, bi,bj,bk,oi,oj,ok);
				  //CpxNumTns& Wblk = W.block(bi,bj,bk);
				  //Wblk(oi,oj,ok) += wpdata(tmpx,tmpy,tmpz)  * bb * sqrt(1.0-ss*ss);
				  O(xcur,ycur,zcur) += wpdata(tmpx,tmpy,tmpz)  * bb * sqrt(1.0-ss*ss);
				}
		  } //ycur
		  
		wcnt++;
	 }
  }//end of face
  //face 5.-z,-x,-y
  for(int g=nd-1; g>=0; g--) {
	 for(int f=nd-1; f>=0; f--) {
		
		  double zs = -R3;		  double ze = -R3/4+(W3/2)/4;
		  double xs = -R1 + (2*f-1)*W1/2;		double xe = -R1 + (2*f+3)*W1/2;
		  double ys = -R2 + (2*g-1)*W2/2;		double ye = -R2 + (2*g+3)*W2/2;
		  int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		  double thts, thtm, thte; //y to x
		  if(f==0) {
			 thts = atan2(-1.0, 1.0-1.0/nd);			 thtm = atan2(-1.0+1.0/nd, 1.0);			 thte = atan2(-1.0+3.0/nd, 1.0);
		  } else if(f==nd-1) {
			 thts = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);			 thtm = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);			 thte = atan2(1.0, 1.0-1.0/nd);
		  } else {
			 thts = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);			 thtm = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);			 thte = atan2(-1.0+(2.0*f+3.0)/nd, 1.0);
		  }
		  double phis, phim, phie; //z to x
		  if(g==0) {
			 phis = atan2(-1.0, 1.0-1.0/nd);			 phim = atan2(-1.0+1.0/nd, 1.0);			 phie = atan2(-1.0+3.0/nd, 1.0);
		  } else if(g==nd-1) {
			 phis = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);			 phim = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);			 phie = atan2(1.0, 1.0-1.0/nd);
		  } else {
			 phis = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);			 phim = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);			 phie = atan2(-1.0+(2.0*g+3.0)/nd, 1.0);
		  }
		  int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2;
		  double R13 = R1/R3;		  double R23 = R2/R3;	  //double R13 = double(F1)/double(F3);		  double R23 = double(F2)/double(F3);
		  CpxNumTns tpdata(xn,yn,zn);
		  //CpxNumTns& Cblk = C.block(s,wcnt);
		  //tpdata = Cblk;
		  ///tpdata = csc[wcnt];
		  tpdata = C.block(s,wcnt);
		  //fft
		  fftwnd_plan p = NULL;
		  map<inttriple, fftwnd_plan>::iterator mit = planmap.find( inttriple(xn, intpair(yn,zn)) );
		  if(mit!=planmap.end()) {			 p = (*mit).second;
		  } else {
			 p = fftw3d_create_plan(zn, yn, xn, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
			 planmap[ inttriple(xn, intpair(yn,zn)) ] = p;
		  }
		  fftwnd_one(p, (fftw_complex*)tpdata.data(), NULL);		  //cerr<<"wedge s"<<endl;
		  double sqrtprod = sqrt(double(xn*yn*zn));
		  for(int i=0; i<xn; i++)			 for(int j=0; j<yn; j++)				for(int k=0; k<zn; k++)				  tpdata(i,j,k) /= sqrtprod;
		  CpxOffTns wpdata(xn,yn,zn);
		  fdct3d_fftshift(xn,yn,zn,tpdata,wpdata);
		  
		  for(int zcur=(int)ceil(zs); zcur<ze; zcur++) {
			 int xfm = (int)ceil( max(-R1, R13*(-zcur)*tan(thts)) );
			 int xto = (int)floor( min(R1, R13*(-zcur)*tan(thte)) );
			 int yfm = (int)ceil( max(-R2, R23*(-zcur)*tan(phis)) );
			 int yto = (int)floor( min(R2, R23*(-zcur)*tan(phie)) );
			 for(int xcur=xfm; xcur<=xto; xcur++)
				for(int ycur=yfm; ycur<=yto; ycur++) {
				  int tmpx = xcur%xn;				  if(tmpx<-xh) tmpx+=xn;				  if(tmpx>=-xh+xn) tmpx-=xn;
				  int tmpy = ycur%yn;				  if(tmpy<-yh) tmpy+=yn;				  if(tmpy>=-yh+yn) tmpy-=yn;
				  int tmpz = zcur%zn;				  if(tmpz<-zh) tmpz+=zn;				  if(tmpz>=-zh+zn) tmpz-=zn;
				  
				  double thtcur = atan2(xcur/R1, (-zcur)/R3);
				  double phicur = atan2(ycur/R2, (-zcur)/R3);
				  double glbpou;					 fdct3d_globalpou(thtcur, phicur, M_PI/4-atan2(1.0-1.0/nd, 1.0), glbpou);
				  double wtht;
				  if(thtcur<thtm) {
					 if(f==0)						wtht = 1;
					 else {						double l,r;						fdct3d_window( (thtcur-thts)/(thtm-thts), l, r);						wtht = l;					 }
				  } else {
					 if(f==nd-1)						wtht = 1;
					 else {						double l,r;						  fdct3d_window( (thtcur-thtm)/(thte-thtm), l, r);						wtht = r;					 }
				  }
				  double wphi;
				  if(phicur<phim) {
					 if(g==0)						wphi = 1;
					 else {						double l,r;						  fdct3d_window( (phicur-phis)/(phim-phis), l, r);						wphi = l;					 }
				  } else {
					 if(g==nd-1)						wphi = 1;
					 else {						double l,r;						  fdct3d_window( (phicur-phim)/(phie-phim), l, r);						wphi = r;					 }
				  }
				  double pou = glbpou * wtht * wphi;
				  wpdata(tmpx, tmpy, tmpz) *= pou;
				  				  
				  double ss = sma1(xcur)*sma2(ycur)*sma3(zcur);				  double bb = big1(xcur)*big2(ycur)*big3(zcur);
				  //int bi,bj,bk;			 int oi,oj,ok;			 fdct3d_position_aux(N1,N2,N3,b, xcur,ycur,zcur, bi,bj,bk,oi,oj,ok);
				  //CpxNumTns& Wblk = W.block(bi,bj,bk);
				  //Wblk(oi,oj,ok) += wpdata(tmpx,tmpy,tmpz)  * bb * sqrt(1.0-ss*ss);
				  O(xcur,ycur,zcur) += wpdata(tmpx,tmpy,tmpz)  * bb * sqrt(1.0-ss*ss);
				}
		  }//zcur

		wcnt++;
	 }
  }//end of face
  assert(wcnt==nd*nd*nf);
  
  //remove plans
  for(map<inttriple, fftwnd_plan>::iterator mit=planmap.begin(); mit!=planmap.end(); mit++) {
	 fftwnd_plan p = (*mit).second;
	 fftwnd_destroy_plan(p);
  }
  
  return 0;
}

int fdct3d_inverse_wavelet(double L1, double L2, double L3, int s, CpxNumTns& W, CpxOffTns& O)
{
  ///vector<CpxNumTns>& csc = C[s];
  
  L1 = L1/2;  L2 = L2/2;  L3 = L3/2;
  int S1, S2, S3;	 int F1, F2, F3;	 double R1, R2, R3;	 fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);
  DblOffVec big1(S1);  fdct3d_lowpass(L1, big1);
  DblOffVec big2(S2);  fdct3d_lowpass(L2, big2);
  DblOffVec big3(S3);  fdct3d_lowpass(L3, big3);
  
  int N1 = O.m();  int N2 = O.n();  int N3 = O.p();
  CpxNumTns T(W);
  fftwnd_plan p = fftw3d_create_plan(N3,N2,N1, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fftwnd_one(p, (fftw_complex*)T.data(), NULL);
  fftwnd_destroy_plan(p);
  double sqrtprod = sqrt(double(N1*N2*N3));
  for(int i=0; i<N1; i++)	 for(int j=0; j<N2; j++)		for(int k=0; k<N3; k++)		  T(i,j,k) /= sqrtprod;
  
  fdct3d_fftshift(N1,N2,N3,T,O);
  
  for(int i=-S1/2; i<-S1/2+S1; i++)
	 for(int j=-S2/2; j<-S2/2+S2; j++)
		for(int k=-S3/2; k<-S3/2+S3; k++) {
		  double pou = big1(i)*big2(j)*big3(k);
		  O(i,j,k) = O(i,j,k) * sqrt(1-pou*pou);
		}
  
  return 0;
}

int fdct3d_inverse_center(double L1, double L2, double L3, int s, CpxCrvletOcr& C, CpxOffTns& O)
{
  //vector<CpxNumTns>& csc = C[s];
  
  int S1, S2, S3;	 int F1, F2, F3;	 double R1, R2, R3;	 fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);
  DblOffVec big1(S1);  fdct3d_lowpass(L1, big1);
  DblOffVec big2(S2);  fdct3d_lowpass(L2, big2);
  DblOffVec big3(S3);  fdct3d_lowpass(L3, big3);

  CpxNumTns T(S1,S2,S3);
  ///T = csc[0];
  T = C.block(s,0);
  fftwnd_plan p = fftw3d_create_plan(S3,S2,S1, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fftwnd_one(p, (fftw_complex*)T.data(), NULL);
  fftwnd_destroy_plan(p);
  double sqrtprod = sqrt(double(S1*S2*S3));
  for(int i=0; i<S1; i++)	 for(int j=0; j<S2; j++)		for(int k=0; k<S3; k++)		  T(i,j,k) /= sqrtprod;
  
  CpxOffTns A(S1,S2,S3);
  fdct3d_fftshift(S1,S2,S3,T,A);
  for(int i=-S1/2; i<-S1/2+S1; i++)
	 for(int j=-S2/2; j<-S2/2+S2; j++)
		for(int k=-S3/2; k<-S3/2+S3; k++) {
		  O(i,j,k) += A(i,j,k) * (big1(i)*big2(j)*big3(k));
		}
 
  return 0;
}


int fdct3d_inverse(int N1, int N2, int N3, int nbscales, int nbdstz_coarse,
					 CpxCrvletOcr& C, CpxNumTns& W,				CpxNumTns& X)
{
  time_t tm0, tm1;
  tm0 = time(NULL); 
  
  X.resize(0,0,0);
  CpxOffTns O(N1,N2,N3);
  
  int L = nbscales;
  {
	 int s = L-1;
	 double L1 = 4.0*N1/3.0 / pow2(L-1-s);	 double L2 = 4.0*N2/3.0 / pow2(L-1-s);	 double L3 = 4.0*N3/3.0 / pow2(L-1-s);
	 fdct3d_inverse_wavelet(L1,L2,L3,s, W,O);
  }//tm1 = time(NULL);  cout<<"wavelet "<<difftime(tm1,tm0)<<endl;  tm0 = tm1;
  {
	 int s = 0;
	 double L1 = 4.0*N1/3.0 / pow2(L-1-s);	 double L2 = 4.0*N2/3.0 / pow2(L-1-s);	 double L3 = 4.0*N3/3.0 / pow2(L-1-s);
	 fdct3d_inverse_center(L1,L2,L3,s, C, O);
  }//tm1 = time(NULL);  cout<<"center "<<difftime(tm1,tm0)<<endl;  tm0 = tm1;
  for(int s=1; s<L-1; s++) {
	 double L1 = 4.0*N1/3.0 / pow2(L-1-s);	 double L2 = 4.0*N2/3.0 / pow2(L-1-s);	 double L3 = 4.0*N3/3.0 / pow2(L-1-s);
	 int nd = nbdstz_coarse * pow2(s/2);
	 fdct3d_inverse_angles(L1,L2,L3,s, nd, C, O);
  }//tm1 = time(NULL);  cout<<"angles "<<difftime(tm1,tm0)<<endl;  tm0 = tm1;
  
  //fft
  CpxNumTns& T = X; T.resize(N1,N2,N3);
  fdct3d_ifftshift(N1,N2,N3, O,T);
  fftwnd_plan p = fftw3d_create_plan(N3, N2, N1, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fftwnd_one(p, (fftw_complex*)T.data(), NULL);
  fftwnd_destroy_plan(p);
  double sqrtprod = sqrt(double(N1*N2*N3)); //scale
  for(int i=0; i<N1; i++)	 for(int j=0; j<N2; j++)		for(int k=0; k<N3; k++)		  T(i,j,k) /= sqrtprod;
  //tm1 = time(NULL);  cout<<"ifft+shift "<<difftime(tm1,tm0)<<endl;  tm0 = tm1;

  return 0;
}
