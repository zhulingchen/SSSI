/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#include "commoninc.hpp"
#include "cpxnumtnsblkd.hpp"
#include "cpxcrvletprtd.hpp"
#include "offvec.hpp"
#include "offmat.hpp"
#include "offtns.hpp"
#include "numvec.hpp"
#include "nummat.hpp"
#include "numtns.hpp"

//----------------------------------------
//forward transform, x->(c,w)
int fdct3d_forward(int N1, int N2, int N3, int nbscales, int nbdstz_coarse,
				   CpxNumTnsBlkd& x,
				   CpxCrvletPrtd& c, CpxNumTnsBlkd& w);
//this function performs the forward curvelet transform
//INPUTS:
//  N1,N2,N3 -- the size of the input data
//  nbscales -- the total number of scales for subband decomposition
//  nbdstz_coarse -- the number of discretizations in each direction on each face in the 2nd coarest scale
//  x -- N1 by N2 by N3 tensor 
//OUTPUTS:
//  c -- the curvelet coefficients data structure coeffient at scale s, wedge w and indices (i,j,k) is accessed by c[s][w](i,j,k).
//  w -- the wavelet coeffients at the finest level, N1 by N2 by N3 tensor

//inverse transform (c,w)->x
int fdct3d_inverse(int N1, int N2, int N3, int nbscales, int nbdstz_coarse,
				   CpxCrvletPrtd& c, CpxNumTnsBlkd& w,
				   CpxNumTnsBlkd& x);
//this function performs the inverse curvelet transform
//INPUTS:
//  N1,N2,N3 -- the size of the input data
//  nbscales -- the total number of scales for subband decomposition
//  nbdstz_coarse -- the number of discretizations in each direction on each face in the 2nd coarest scale
//  c -- the curvelet coefficients data structure coeffient at scale s, wedge w and indices (i,j,k) is accessed by c[s][w](i,j,k). 
//  w -- the wavelet coeffients at the finest level, N1 by N2 by N3 tensor
//OUTPUTS:
//  x -- N1 by N2 by N3 tensor 


//for the curvelet transform, extract the freqency center and grid size for each curvelet
int fdct3d_param(  int N1, int N2, int N3, int nbscales, int nbdstz_coarse, 
				   vector< vector<double> >& fxs, vector< vector<double> >& fys, vector< vector<double> >& fzs,
				   vector< vector<int   > >& nxs, vector< vector<int   > >& nys, vector< vector<int   > >& nzs);
//this function obtains auxiliary information about curvelet transform 
//INPUTS:
//  N1,N2,N3 -- the size of the input data
//  nbscales -- the total number of scales for subband decomposition
//  nbdstz_coarse -- the number of discretizations in each direction on each face in the 2nd coarest scale
//OUTPUTS:
//  fx,fy,fz -- for scale s and wege w, fx[s][w], fy[s][w] and fz[s][w] give the coordinate of the frequency center of curvelet of scale s and wedge w
//  nx,ny,nz -- for scale s and wege w, nx[s][w], ny[s][w] and nz[s][w] give the dimensions of the curvelet matrix of scale s and wedge w


//----------------------------------------
//auxilary functions
int fdct3d_partition_cpxnumtnsblkd_z(int N1,int N2,int N3,int b, 
									 BolNumTns& exists, IntNumTns& owners);
int fdct3d_partition_cpxnumtnsblkd_y(int N1,int N2,int N3,int b, 
									 BolNumTns& exists, IntNumTns& owners);
int fdct3d_partition_cpxnumtnsblkd_x(int N1,int N2,int N3,int b, 
									 BolNumTns& exists, IntNumTns& owners);

int fdct3d_partition_cpxcrvletprtd(int N1,int N2,int N3, int nbscales,int nbdstz_coarse,
								   vector< vector<bool> >& exists, vector< vector<int> >& owners);

int fdct3d_dependency(int N1,int N2,int N3, int b, int nbscales, int nbdstz_coarse,	vector< vector<int> >& crvowners, 
					  BolNumTns& tnsexists);

int fdct3d_fft( CpxNumTnsBlkd& W);

int fdct3d_ifft(CpxNumTnsBlkd& W);

