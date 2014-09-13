/*
  Copyright (C) 2004 Caltech
  Written by Lexing Ying
*/

#ifndef _FDCT_USFFT_INC_HPP_
#define _FDCT_USFFT_INC_HPP_

#include <iostream>
#include <fstream>
#include <sstream>

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
#include <string.h>
//FFT stuff
#include "fftw.h"

#define FDCT_USFFT_NS_BEGIN_NAMESPACE namespace fdct_usfft_ns {
#define FDCT_USFFT_NS_END_NAMESPACE }

FDCT_USFFT_NS_BEGIN_NAMESPACE

typedef std::complex<double> cpx;

//AUX functions
inline int pow2(int l) { assert(l>=0); return (1<<l); }

FDCT_USFFT_NS_END_NAMESPACE

#endif
