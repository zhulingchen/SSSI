/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#ifndef _FDCT3DINC_HPP_
#define _FDCT3DINC_HPP_

//STL stuff
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <cfloat>
#include <cassert>
#include <cmath>
#include <string>
#include <complex>

#include <vector>
#include <set>
#include <map>
#include <deque>
#include <queue>
#include <utility>
#include <algorithm>
using namespace std;

//FFT stuff
#include "fftw.h"

//typedef double double;
typedef complex<double> cpx;

//AUX functions
inline int pow2(int l) { assert(l>=0); return (1<<l); }

#endif
