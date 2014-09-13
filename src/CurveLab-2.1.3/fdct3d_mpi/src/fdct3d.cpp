/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#include "fdct3d.hpp"
#include "fdct3dinline.hpp"

//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
int fdct3d_partition_cpxnumtnsblkd_z(int N1,int N2,int N3,int b, 
									 BolNumTns& exists, IntNumTns& owners)
{
  int e = N1/b;  int f = N2/b;  int g = N3/b;
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  iA(g%mpisize==0);
  int bpp = g/mpisize; //block per processor
  for(int k=0; k<g; k++) {
	 int pi = k/bpp; //owner processor
	 for(int i=0; i<e; i++)		for(int j=0; j<f; j++) {
		owners(i,j,k) = pi;
		exists(i,j,k) = (pi==mpirank);
	 }
  }
  return 0;
}

int fdct3d_partition_cpxnumtnsblkd_y(int N1,int N2,int N3,int b, 
									 BolNumTns& exists, IntNumTns& owners)
{
  int e = N1/b;  int f = N2/b;  int g = N3/b;
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  iA(f%mpisize==0);
  int bpp = f/mpisize;
  for(int j=0; j<f; j++) {
	 int pi = j/bpp; //owner processor
	 for(int i=0; i<e; i++)		for(int k=0; k<g; k++) {
		owners(i,j,k) = pi;
		exists(i,j,k) = (pi==mpirank);
	 }
  }
  return 0;
}

int fdct3d_partition_cpxnumtnsblkd_x(int N1,int N2,int N3,int b, 
									 BolNumTns& exists, IntNumTns& owners)
{
  int e = N1/b;  int f = N2/b;  int g = N3/b;
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  iA(e%mpisize==0);
  int bpp = e/mpisize;
  for(int i=0; i<e; i++) {
	 int pi = i/bpp; //owner processor
	 for(int j=0; j<f; j++)		for(int k=0; k<g; k++) {
		owners(i,j,k) = pi;
		exists(i,j,k) = (pi==mpirank);
	 }
  }
  return 0;
}

//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
int fdct3d_partition_cpxcrvletprtd(int N1,int N2,int N3, int nbscales,int nbdstz_coarse,
								   vector< vector<bool> >& exists, vector< vector<int> >& owners)
{
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  
  exists.resize(nbscales-1);
  owners.resize(nbscales-1);
  {
	 int s = 0;	 vector<bool>& es = exists[s];	 vector<int>& os = owners[s];
	 es.resize(1);	 os.resize(1);
	 int pi = 0;
	 os[0] = pi;
	 es[0] = (pi==mpirank);
  }
  int nf = 6;
  
  for(int s=1; s<nbscales-1; s++) {
	 int nd = nbdstz_coarse * pow2(s/2);
	 vector<bool>& es = exists[s];	 vector<int>& os = owners[s];
	 es.resize(nd*nd * nf);
	 os.resize(nd*nd * nf);
	 
	 int face;
	 int wcnt = 0;
	 int wpp = (nd*nd*nf)/mpisize;	 if(wpp*mpisize<nd*nd*nf) wpp++;
	 //face 0, x y z
	 face = 0;
	 for(int h=0; h<nd; h++) {
		for(int g=0; g<nd; g++) {
		  int pi = wcnt/wpp;
		  os[wcnt] = pi;
		  es[wcnt] = (pi==mpirank);
		  wcnt++;
		}
	 }
	 //face 1 y z x
	 face = 1;
	 for(int f=0; f<nd; f++) {
		for(int h=0; h<nd; h++) {
		  int pi = wcnt/wpp;
		  os[wcnt] = pi;
		  es[wcnt] = (pi==mpirank);
		  wcnt++;
		}
	 }
	 //face 2 z x y
	 face = 2;
	 for(int g=0; g<nd; g++) {
		for(int f=0; f<nd; f++) {
		  int pi = wcnt/wpp;
		  os[wcnt] = pi;
		  es[wcnt] = (pi==mpirank);
		  wcnt++;
		}
	 }
	 //face 3 -x -y -z
	 face = 3;
	 for(int h=nd-1; h>=0; h--) {
		for(int g=nd-1; g>=0; g--) {
		  int pi = wcnt/wpp;
		  os[wcnt] = pi;
		  es[wcnt] = (pi==mpirank);
		  wcnt++;
		}
	 }
	 //face 4 -y -z -x
	 face = 4;
	 for(int f=nd-1; f>=0; f--) {
		for(int h=nd-1; h>=0; h--) {
		  int pi = wcnt/wpp;
		  os[wcnt] = pi;
		  es[wcnt] = (pi==mpirank);
		  wcnt++;
		}
	 }
	 //face 5 -z -x -y
	 face = 5;
	 for(int g=nd-1; g>=0; g--) {
		for(int f=nd-1; f>=0; f--) {
		int pi = wcnt/wpp;
		  os[wcnt] = pi;
		  es[wcnt] = (pi==mpirank);
		  wcnt++;
		}
	 }
	 iA(wcnt==nf*nd*nd);
  } //scale
  
  return 0;
}

//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
int fdct3d_dependency_center(int N1,int N2,int N3, int b, double L1, double L2, double L3, vector<int>& crvowners,
								  BolNumTns& tnsexists)
{
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  
  int S1, S2, S3;	 int F1, F2, F3;	 double R1, R2, R3;	 fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);
  if(crvowners[0]==mpirank) {
	 for(int xcur=-F1; xcur<=F1; xcur++)
		for(int ycur=-F2; ycur<=F2; ycur++)
		  for(int zcur=-F3; zcur<=F3; zcur++) {
			 int bi,bj,bk;			 int oi,oj,ok;			 fdct3d_position_aux(N1,N2,N3,b,xcur,ycur,zcur,bi,bj,bk,oi,oj,ok);
			 tnsexists(bi,bj,bk) = true;
		  }
  }
  return 0;
}

int fdct3d_dependency_angles(int N1,int N2,int N3, int b, double L1, double L2, double L3, int nd, vector<int>& crvowners,
								  BolNumTns& tnsexists)
{
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  
  int nf = 6;
  int wcnt = 0;
  int S1, S2, S3;	 int F1, F2, F3;	 double R1, R2, R3;	 fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);
  double W1 = L1/nd;  double W2 = L2/nd;  double W3 = L3/nd;
  //face 0: x,y,z
  for(int h=0; h<nd; h++) {
	 for(int g=0; g<nd; g++) {
		if(crvowners[wcnt]==mpirank) {
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
		  //int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2;
		  double R21 = R2/R1;		  double R31 = R3/R1;
		  for(int xcur=(int)ceil(xs); xcur<xe; xcur++) {
			 int yfm = (int)ceil( max(-R2, R21*xcur*tan(thts)) );
			 int yto = (int)floor( min(R2, R21*xcur*tan(thte)) );
			 int zfm = (int)ceil( max(-R3, R31*xcur*tan(phis)) );
			 int zto = (int)floor( min(R3, R31*xcur*tan(phie)) );
			 for(int ycur=yfm; ycur<=yto; ycur++)
				for(int zcur=zfm; zcur<=zto; zcur++) {
				  int bi,bj,bk;			 int oi,oj,ok;			 fdct3d_position_aux(N1,N2,N3,b,xcur,ycur,zcur,bi,bj,bk,oi,oj,ok);
				  tnsexists(bi,bj,bk) = true;
				}
		  } //xcur
		} //if
		wcnt++;
	 }
  } //end of face
  //face 1: y,z,x
  for(int f=0; f<nd; f++) {
	 for(int h=0; h<nd; h++) {
		if(crvowners[wcnt]==mpirank) {
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
		  //int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2;
		  double R32 = R3/R2;		  double R12 = R1/R2;
		  for(int ycur=(int)ceil(ys); ycur<ye; ycur++) {
			 int zfm = (int)ceil( max(-R3, R32*ycur*tan(thts)) );
			 int zto = (int)floor( min(R3, R32*ycur*tan(thte)) );
			 int xfm = (int)ceil( max(-R1, R12*ycur*tan(phis)) );
			 int xto = (int)floor( min(R1, R12*ycur*tan(phie)) );
			 for(int zcur=zfm; zcur<=zto; zcur++)
				for(int xcur=xfm; xcur<=xto; xcur++) {
				  int bi,bj,bk;			 int oi,oj,ok;			 fdct3d_position_aux(N1,N2,N3,b,xcur,ycur,zcur,bi,bj,bk,oi,oj,ok);
				  tnsexists(bi,bj,bk) = true;
				}
		  } //ycur
		}//if
		wcnt++;
	 }
  } //end of face
  //face 2: z,x,y
  for(int g=0; g<nd; g++) {
	 for(int f=0; f<nd; f++) {
		if(crvowners[wcnt]==mpirank) {
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
		  //int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2;
		  double R13 = double(F1)/double(F3);		  double R23 = double(F2)/double(F3);
		  for(int zcur=(int)ceil(zs); zcur<ze; zcur++) {
			 int xfm = (int)ceil( max(-R1, R13*zcur*tan(thts)) );
			 int xto = (int)floor( min(R1, R13*zcur*tan(thte)) );
			 int yfm = (int)ceil( max(-R2, R23*zcur*tan(phis)) );
			 int yto = (int)floor( min(R2, R23*zcur*tan(phie)) );
			 for(int xcur=xfm; xcur<=xto; xcur++)
				for(int ycur=yfm; ycur<=yto; ycur++) {
				  int bi,bj,bk;			 int oi,oj,ok;			 fdct3d_position_aux(N1,N2,N3,b,xcur,ycur,zcur,bi,bj,bk,oi,oj,ok);
				  tnsexists(bi,bj,bk) = true;
				}
		  }//zcur
		}//if
		wcnt++;
	 }
  } //end of face
  //face 3: -x,-y,-z
  for(int h=nd-1; h>=0; h--) {
	 for(int g=nd-1; g>=0; g--) {
		if(crvowners[wcnt]==mpirank) {
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
		  //int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2;
		  double R21 = R2/R1;		  double R31 = R3/R1;
		  for(int xcur=(int)ceil(xs); xcur<xe; xcur++) {
			 int yfm = (int)ceil( max(-R2, R21*(-xcur)*tan(thts)) );
			 int yto = (int)floor( min(R2, R21*(-xcur)*tan(thte)) );
			 int zfm = (int)ceil( max(-R3, R31*(-xcur)*tan(phis)) );
			 int zto = (int)floor( min(R3, R31*(-xcur)*tan(phie)) );
			 for(int ycur=yfm; ycur<=yto; ycur++)
				for(int zcur=zfm; zcur<=zto; zcur++) {
				  int bi,bj,bk;			 int oi,oj,ok;			 fdct3d_position_aux(N1,N2,N3,b,xcur,ycur,zcur,bi,bj,bk,oi,oj,ok);
				  tnsexists(bi,bj,bk) = true;
				}
		  } //xcur
		} //if
		wcnt++;
	 }
  } //end of face
  //face 4: -y,-z,-x
  for(int f=nd-1; f>=0; f--) {
	 for(int h=nd-1; h>=0; h--) {
		if(crvowners[wcnt]==mpirank) {
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
		  //int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2;
		  double R32 = double(F3)/double(F2);		  double R12 = double(F1)/double(F2);
		  for(int ycur=(int)ceil(ys); ycur<ye; ycur++) {
			 int zfm = (int)ceil( max(-R3, R32*(-ycur)*tan(thts)) );
			 int zto = (int)floor( min(R3, R32*(-ycur)*tan(thte)) );
			 int xfm = (int)ceil( max(-R1, R12*(-ycur)*tan(phis)) );
			 int xto = (int)floor( min(R1, R12*(-ycur)*tan(phie)) );
			 for(int zcur=zfm; zcur<=zto; zcur++)
				for(int xcur=xfm; xcur<=xto; xcur++) {
				  int bi,bj,bk;			 int oi,oj,ok;			 fdct3d_position_aux(N1,N2,N3,b,xcur,ycur,zcur,bi,bj,bk,oi,oj,ok);
				  tnsexists(bi,bj,bk) = true;
				}
		  } //ycur
		}//if
		wcnt++;
	 }
  }//end of face
  //face 5: -z,-x,-y
  for(int g=nd-1; g>=0; g--) {
	 for(int f=nd-1; f>=0; f--) {
		if(crvowners[wcnt]==mpirank) {
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
		  //int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2;
		  double R13 = double(F1)/double(F3);		  double R23 = double(F2)/double(F3);
		  for(int zcur=(int)ceil(zs); zcur<ze; zcur++) {
			 int xfm = (int)ceil( max(-R1, R13*(-zcur)*tan(thts)) );
			 int xto = (int)floor( min(R1, R13*(-zcur)*tan(thte)) );
			 int yfm = (int)ceil( max(-R2, R23*(-zcur)*tan(phis)) );
			 int yto = (int)floor( min(R2, R23*(-zcur)*tan(phie)) );
			 for(int xcur=xfm; xcur<=xto; xcur++)
				for(int ycur=yfm; ycur<=yto; ycur++) {
				  int bi,bj,bk;			 int oi,oj,ok;			 fdct3d_position_aux(N1,N2,N3,b,xcur,ycur,zcur,bi,bj,bk,oi,oj,ok);
				  tnsexists(bi,bj,bk) = true;
				}
		  }//zcur
		}//if
		wcnt++;
	 }
  }//end of face
  iA(wcnt==nd*nd*nf);
  
  return 0;
}
int fdct3d_dependency(int N1,int N2,int N3, int b, int nbscales, int nbdstz_coarse,	vector< vector<int> >& crvowners, 
						 BolNumTns& tnsexists)
{
  //int e = N1/b;  int f = N2/b;  int g = N3/b;
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  
  setvalue(tnsexists, false);
  int L = nbscales;
  {
	 int s = 0;
	 double L1 = 2.0*N1/3.0 / pow2(L-2-s);	 double L2 = 2.0*N2/3.0 / pow2(L-2-s);	 double L3 = 2.0*N3/3.0 / pow2(L-2-s);
	 fdct3d_dependency_center(N1,N2,N3,b, L1,L2,L3, crvowners[s], tnsexists);
  }
  for(int s=1; s<nbscales-1; s++) {
	 double L1 = 2.0*N1/3.0 / pow2(L-2-s);	 double L2 = 2.0*N2/3.0 / pow2(L-2-s);	 double L3 = 2.0*N3/3.0 / pow2(L-2-s);
	 int nd = nbdstz_coarse * pow2(s/2);
	 fdct3d_dependency_angles(N1,N2,N3,b, L1,L2,L3, nd, crvowners[s], tnsexists);
  }
  return 0;
}


//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
int fdct3d_fft(CpxNumTnsBlkd& W)
{
  //LEXING: entering with z slices, exiting with y slices, and the frequency is centered
  time_t tm0, tm1;  tm0 = time(NULL); 
  
  int m = W.m();  int n = W.n();  int p = W.p();
  int N1 = m;  int N2 = n;  int N3 = p;
  int b = W.b();
  int e = W.e();  int f = W.f();  int g = W.g();
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  IntNumTns& Wowners = W.owners();
  double sqrtprod = sqrt(double(N1*N2*N3));
  
  //scale w's value by a complex exponent
  for(int i=0; i<e; i++)	 for(int j=0; j<f; j++)		for(int k=0; k<g; k++) {
	 if(Wowners(i,j,k)==mpirank) { //i am the owner of this block
		CpxNumTns& Wblk = W.block(i,j,k);
		int istt = i*b;		int jstt = j*b;		int kstt = k*b;
		for(int ioff=0; ioff<b; ioff++)		  for(int joff=0; joff<b; joff++)			 for(int koff=0; koff<b; koff++) {
		  cpx coef = polar( 1.0, M_PI*(ioff+istt + joff+jstt + koff+kstt) );
		  Wblk(ioff, joff, koff) *= coef/sqrtprod;
		}
	 }
  }
  //copy w's z slices into (x,y) cpxnummat and perform fft on each of them, then put back
  fftwnd_plan pxyf = fftw2d_create_plan(N2, N1, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
  for(int k=0; k<g; k++) {
	 int kstt = k*b;
	 if(Wowners(0,0,k)==mpirank) { //i am the owner of this z slice
		//for each silce, store, fft and put back
		for(int koff=0; koff<b; koff++) { 
		  CpxNumMat tmp(N1, N2);
		  for(int i=0; i<e; i++)		  for(int j=0; j<f; j++) {
			 iA( Wowners(i,j,k)==mpirank); //make sure i am the owner
			 CpxNumTns& Wblk = W.block(i,j,k);
			 int istt = i*b;		int jstt = j*b;
			 for(int ioff=0; ioff<b; ioff++)			 for(int joff=0; joff<b; joff++)
				tmp(ioff+istt, joff+jstt) = Wblk(ioff,joff,koff);
		  }
		  //fft it
		  fftwnd_one(pxyf, (fftw_complex*)(tmp.data()), NULL);
		  //put back
		  for(int i=0; i<e; i++)			 for(int j=0; j<f; j++) {
			 iA( Wowners(i,j,k)==mpirank);
			 CpxNumTns& Wblk = W.block(i,j,k);
			 int istt = i*b;			 int jstt = j*b;
			 for(int ioff=0; ioff<b; ioff++)			 for(int joff=0; joff<b; joff++)
				Wblk(ioff,joff,koff) = tmp(ioff+istt, joff+jstt);
		  }
		} //koff
	 }
  }
  fftwnd_destroy_plan(pxyf);  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "  fft compute1 %f\n", difftime(tm1,tm0)) );  tm0 = tm1;

  BolNumTns newtnsexists(e,f,g);
  IntNumTns newtnsowners(e,f,g);
  //scatter w to contain y slices
  fdct3d_partition_cpxnumtnsblkd_y(m,n,p,b,newtnsexists,newtnsowners);
  iC( W.scatter(newtnsexists) );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "  fft scatter %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  //shift w's owner to y slices
  iC( W.shift(newtnsowners) );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "  fft shift %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  //discard w's nonowners
  iC( W.discard() );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "  fft discard %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  
  //copy w's data into (z) cpxnumvec and perform fft on each of them, then put back
  fftw_plan pzf = fftw_create_plan(N3, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
  for(int i=0; i<e; i++)  for(int j=0; j<f; j++) {
	 int istt = i*b;	 int jstt = j*b;
	 if(Wowners(i,j,0)==mpirank) {
		//for each line, store, fft and put back
		for(int ioff=0; ioff<b; ioff++)		  for(int joff=0; joff<b; joff++) {
		  CpxNumVec tmp(N3);
		  for(int k=0; k<g; k++) {
			 iA( Wowners(i,j,k)==mpirank );
			 CpxNumTns& Wblk = W.block(i,j,k);
			 int kstt = k*b;
			 for(int koff=0; koff<b; koff++)
				tmp(koff+kstt) = Wblk(ioff, joff, koff);
		  }
		  //fft it
		  fftw_one(pzf, (fftw_complex*)(tmp.data()), NULL);
		  //put back
		  for(int k=0; k<g; k++) {
			 iA( Wowners(i,j,k)==mpirank );
			 CpxNumTns& Wblk = W.block(i,j,k);
			 int kstt = k*b;
			 for(int koff=0; koff<b; koff++)
				Wblk(ioff,joff,koff) = tmp(koff+kstt);
		  }
		}
	 }
  }
  fftw_destroy_plan(pzf);  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "  fft compute2 %f\n", difftime(tm1,tm0)) );  tm0 = tm1;

  return 0;
}

//-----------------------------------------------------------------
int fdct3d_ifft(CpxNumTnsBlkd& W)
{
  //LEXING: entering with y slices, exiting with z slices, and the frequency is centered
  time_t tm0, tm1;  tm0 = time(NULL); 
  
  int m = W.m();  int n = W.n();  int p = W.p();
  int N1 = m;  int N2 = n;  int N3 = p;
  int b = W.b();
  int e = W.e();  int f = W.f();  int g = W.g();
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  IntNumTns& Wowners = W.owners();
  double sqrtprod = sqrt(double(N1*N2*N3));
  
  //copy w's data into (z) cpxnumvec and perform ifft on each of them, then put back
  fftw_plan pzb = fftw_create_plan(N3, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
  for(int i=0; i<e; i++)	 for(int j=0; j<f; j++) {
	 int istt = i*b;	 int jstt = j*b;
	 if(Wowners(i,j,0)==mpirank) {
		//for each line, store, ifft and put back
		for(int ioff=0; ioff<b; ioff++)		  for(int joff=0; joff<b; joff++) {
		  CpxNumVec tmp(N3);
		  for(int k=0; k<g; k++) {
			 iA( Wowners(i,j,k)==mpirank );
			 CpxNumTns& Wblk = W.block(i,j,k);
			 int kstt = k*b;
			 for(int koff=0; koff<b; koff++)
				tmp(koff+kstt) = Wblk(ioff, joff, koff);
		  }
		  //fft it
		  fftw_one(pzb, (fftw_complex*)(tmp.data()), NULL);
		  //put back
		  for(int k=0; k<g; k++) {
			 iA( Wowners(i,j,k)==mpirank);
			 CpxNumTns& Wblk = W.block(i,j,k);
			 int kstt = k*b;
			 for(int koff=0; koff<b; koff++)
				Wblk(ioff,joff,koff) = tmp(koff+kstt);
		  }
		}
	 }
  }
  fftw_destroy_plan(pzb);  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "  ifft compute1 %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  
  BolNumTns newtnsexists(e,f,g);
  IntNumTns newtnsowners(e,f,g);
  //scatter w to contain z slices
  fdct3d_partition_cpxnumtnsblkd_z(m,n,p,b,newtnsexists,newtnsowners);
  iC( W.scatter(newtnsexists) );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "  ifft scatter %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  //shift w's onwer to z slices
  iC( W.shift(newtnsowners) );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "  ifft shift %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  //discard w's nonowners
  iC( W.discard() );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "  ifft discard %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  
  //copy w's data into (x,y) cpxnummat and perform ifft on each of them, then put back
  fftwnd_plan pxyb = fftw2d_create_plan(N2, N1, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
  for(int k=0; k<g; k++) {
	 int kstt = k*b;
	 if(Wowners(0,0,k)==mpirank) {
		for(int koff=0; koff<b; koff++) {
		  CpxNumMat tmp(N1, N2);
		  for(int i=0; i<e; i++)		  for(int j=0; j<f; j++) {
			 iA( Wowners(i,j,k)==mpirank); //make sure i am the owner
			 CpxNumTns& Wblk = W.block(i,j,k);
			 int istt = i*b;		int jstt = j*b;
			 for(int ioff=0; ioff<b; ioff++)			 for(int joff=0; joff<b; joff++)
				tmp(ioff+istt, joff+jstt) = Wblk(ioff,joff,koff);
		  }
		  //fft it
		  fftwnd_one(pxyb, (fftw_complex*)(tmp.data()), NULL);
		  //put back
		  for(int i=0; i<e; i++)			 for(int j=0; j<f; j++) {
			 iA( Wowners(i,j,k)==mpirank);
			 CpxNumTns& Wblk = W.block(i,j,k);
			 int istt = i*b;			 int jstt = j*b;
			 for(int ioff=0; ioff<b; ioff++)			 for(int joff=0; joff<b; joff++)
				Wblk(ioff,joff,koff) = tmp(ioff+istt, joff+jstt);
		  }
		}		
	 }
  }
  fftwnd_destroy_plan(pxyb);  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "  ifft compute2 %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  
  //scale w's value by a complex exponent
  for(int i=0; i<e; i++)	 for(int j=0; j<f; j++)		for(int k=0; k<g; k++) {
	 if(Wowners(i,j,k)==mpirank) { //i am the owner of this block
		CpxNumTns& Wblk = W.block(i,j,k);
		int istt = i*b;		int jstt = j*b;		int kstt = k*b;
		for(int ioff=0; ioff<b; ioff++)		  for(int joff=0; joff<b; joff++)			 for(int koff=0; koff<b; koff++) {
		  cpx coef = polar( 1.0, - M_PI*(ioff+istt + joff+jstt + koff+kstt) ); //negative size
		  Wblk(ioff, joff, koff) *= coef/sqrtprod;
		}
	 }
  }
  
  return 0;
}


