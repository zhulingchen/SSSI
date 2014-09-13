/*
  Copyright (C) 2004 Caltech
  Written by Lexing Ying
*/

#ifndef _FDCT_USFFT_HPP_
#define _FDCT_USFFT_HPP_

#include "fdct_usfft_inc.hpp"
#include "numvec.hpp"
#include "nummat.hpp"
#include "offvec.hpp"
#include "offmat.hpp"
using std::vector;
using std::map;
using std::sort;
using std::pair;
using std::max;
using std::min;
using std::abs;
using std::polar;

FDCT_USFFT_NS_BEGIN_NAMESPACE

struct dbl2 {
  double _x, _y;
  dbl2(): _x(0), _y(0) {;}
  dbl2(double x, double y): _x(x), _y(y) {;}
};

int fdct_usfft( int N1, int N2, int nbscales, int nbangles_coarse, int ac, CpxNumMat& x, vector< vector<CpxNumMat> >& c);
//this function performs the forward curvelet transform
//INPUTS:
//  N1,N2 -- the size of the input image
//  nbscales -- the total number of scales for subband decomposition
//  nbangles_coarse -- the number of angles in the 2nd coarest scale
//  ac -- ac==1 use curvelet at the finest level, ac==0, use wavelet at the finest level 
//  x -- N1 by N2 matrix stored in CpxNumMat class
//OUTPUTS:
//  c -- the curvelet coefficients data structure coeffient at scale s, wedge w and indices (i,j) is accessed by c[s][w](i,j). 

int afdct_usfft(int N1, int N2, int nbscales, int nbangles_coarse, int ac, vector< vector<CpxNumMat> >& c, CpxNumMat& x);
//this function performs the adjoint curvelet transform
//INPUTS:
//  N1,N2 -- the size of the input image
//  nbscales -- the total number of scales for subband decomposition
//  nbangles_coarse -- the number of angles in the 2nd coarest scale
//  ac -- ac==1 use curvelet at the finest level, ac==0, use wavelet at the finest level 
//  c -- the curvelet coefficients data structure coeffient at scale s, wedge w and indices (i,j) is accessed by c[s][w](i,j). 
//OUTPUTS:
//  x -- N1 by N2 matrix stored in CpxNumMat class

int ifdct_usfft(int N1, int N2, int nbscales, int nbangles_coarse, int ac, vector< vector<CpxNumMat> >& c, CpxNumMat& x);
//this function performs the inverse curvelet transform
//INPUTS:
//  N1,N2 -- the size of the input image
//  nbscales -- the total number of scales for subband decomposition
//  nbangles_coarse -- the number of angles in the 2nd coarest scale
//  ac -- ac==1 use curvelet at the finest level, ac==0, use wavelet at the finest level 
//  c -- the curvelet coefficients data structure coeffient at scale s, wedge w and indices (i,j) is accessed by c[s][w](i,j). 
//OUTPUTS:
//  x -- N1 by N2 matrix stored in CpxNumMat class

int fdct_usfft_param(int N1, int N2, int nbscales, int nbangles_coarse, int ac,
							vector< vector<dbl2> >& sx, vector< vector<dbl2> >& sy,
							vector< vector<double> >& fx, vector< vector<double> >& fy,
							vector< vector<int> >& nx, vector< vector<int> >& ny);
//this function obtains auxiliary information about curvelet transform 
//INPUTS:
//  N1,N2 -- the size of the input image
//  nbscales -- the total number of scales for subband decomposition
//  nbangles_coarse -- the number of angles in the 2nd coarest scale
//  ac -- ac==1 use curvelet at the finest level, ac==0, use wavelet at the finest level 
//OUTPUTS:
//  sx -- for scale s and wege w, sx[s][w] gives two doubles denoting the difference vector of the spatial centers of two 
//        consecutive curvelets with the same column index (i.e. (i,j) and (i+1,j))
//  sx -- for scale s and wege w, sy[s][w] gives two doubles denoting the difference vector of the spatial centers of two 
//        consecutive curvelets with the same row index (i.e. (i,j) and (i,j+1))
//  fx -- for scale s and wege w, fx[s][w] is the first coordinate of the frequency center of curvelet of scale s and wedge w
//  fy -- for scale s and wege w, fy[s][w] is the second coordinate of the frequency center of curvelet of scale s and wedge w 
//  nx -- for scale s and wege w, nx[s][w] is the first dimension of the curvelet matrix of scale s and wedge w
//  ny -- for scale s and wege w, ny[s][w] is the second dimension of the curvelet matrix of scale s and wedge w 

FDCT_USFFT_NS_END_NAMESPACE

#endif
