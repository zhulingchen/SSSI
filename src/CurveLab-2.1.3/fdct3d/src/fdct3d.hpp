/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/


#ifndef _FDCT3D_HPP_
#define _FDCT3D_HPP_

#include "commoninc.hpp"
#include "numvec.hpp"
#include "nummat.hpp"
#include "numtns.hpp"
#include "offvec.hpp"
#include "offmat.hpp"
#include "offtns.hpp"

int fdct3d_forward(int N1, int N2, int N3, int nbscales, int nbdstz_coarse, int ac, CpxNumTns& x, vector< vector<CpxNumTns> >& c);
//this function performs the forward curvelet transform
//INPUTS:
//  N1,N2,N3 -- the size of the input data
//  nbscales -- the total number of scales for subband decomposition
//  nbdstz_coarse -- the number of discretizations in each direction on each face in the 2nd coarest scale
//  ac -- ac==1 use curvelet at the finest level, ac==0, use wavelet at the finest level 
//  x -- N1 by N2 by N3 tensor stored in CpxNumTns class
//OUTPUTS:
//  c -- the curvelet coefficients data structure coeffient at scale s, wedge w and indices (i,j,k) is accessed by c[s][w](i,j,k).


int fdct3d_inverse(int N1, int N2, int N3, int nbscales, int nbdstz_coarse, int ac, vector< vector<CpxNumTns> >& c, CpxNumTns& x);
//this function performs the inverse curvelet transform
//INPUTS:
//  N1,N2,N3 -- the size of the input data
//  nbscales -- the total number of scales for subband decomposition
//  nbdstz_coarse -- the number of discretizations in each direction on each face in the 2nd coarest scale
//  ac -- ac==1 use curvelet at the finest level, ac==0, use wavelet at the finest level 
//  c -- the curvelet coefficients data structure coeffient at scale s, wedge w and indices (i,j,k) is accessed by c[s][w](i,j,k). 
//OUTPUTS:
//  x -- N1 by N2 by N3 tensor stored in CpxNumTns class


int fdct3d_param(  int N1, int N2, int N3, int nbscales, int nbdstz_coarse, int ac,
						 vector< vector<double> >& fxs, vector< vector<double> >& fys, vector< vector<double> >& fzs,
						 vector< vector<int   > >& nxs, vector< vector<int   > >& nys, vector< vector<int   > >& nzs);
//this function obtains auxiliary information about curvelet transform 
//INPUTS:
//  N1,N2,N3 -- the size of the input data
//  nbscales -- the total number of scales for subband decomposition
//  nbdstz_coarse -- the number of discretizations in each direction on each face in the 2nd coarest scale
//  ac -- ac==1 use curvelet at the finest level, ac==0, use wavelet at the finest level 
//OUTPUTS:
//  fx,fy,fz -- for scale s and wege w, fx[s][w], fy[s][w] and fz[s][w] give the coordinate of the frequency center of curvelet of scale s and wedge w
//  nx,ny,nz -- for scale s and wege w, nx[s][w], ny[s][w] and nz[s][w] give the dimensions of the curvelet matrix of scale s and wedge w


#endif
