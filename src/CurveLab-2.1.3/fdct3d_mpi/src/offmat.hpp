/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#ifndef _OFFMAT_HPP_
#define _OFFMAT_HPP_

#include "commoninc.hpp"

//-------------------------------------------------
template <class F>
class OffMat {
public:
  int _m, _n;
  int _s, _t;
  F* _data;
public:
  OffMat(int m=0, int n=0): _m(m), _n(n), _s(-m/2), _t(-n/2) {
	 if(_m>0 && _n>0) {		_data = new F[_m*_n]; assert( _data!=NULL );	 memset(_data, 0, _m*_n*sizeof(F));	 } else		_data = NULL;
  }
  OffMat(const OffMat& C): _m(C._m), _n(C._n), _s(C._s), _t(C._t) {
	 if(_m>0 && _n>0) {		_data = new F[_m*_n]; assert( _data!=NULL );	 memset(_data, 0, _m*_n*sizeof(F));	 } else 		_data = NULL;
	 if(_m>0 && _n>0) {		memcpy( _data, C._data, _m*_n*sizeof(F) );	 }
  }
  ~OffMat() {
	 if(_m>0 && _n>0) {		delete[] _data; _data = NULL; }
  }
  OffMat& operator=(const OffMat& C) {
	 if(_m>0 && _n>0) {		delete[] _data; _data = NULL; }
	 _m = C._m; _n = C._n;	 _s = C._s; _t = C._t;
	 if(_m>0 && _n>0) {		_data = new F[_m*_n]; assert( _data!=NULL );	 memset(_data, 0, _m*_n*sizeof(F));	 } else		_data = NULL;
	 if(_m>0 && _n>0) {		memcpy( _data, C._data, _m*_n*sizeof(F) );	 }
	 return *this;
  }
  void resize(int m, int n)  {
	 if(_m!=m || _n!=n) {
		if(_m>0 && _n>0) {		delete[] _data; _data = NULL;		}
		_m = m; _n = n; _s = -m/2; _t = -n/2;
		if(_m>0 && _n>0) {		_data = new F[_m*_n]; assert( _data!=NULL );		memset(_data, 0, _m*_n*sizeof(F));		} else		  _data = NULL;
	 }
  }
  const F& operator()(int i, int j) const  { 	 assert( i>=_s && i<_m+_s && j>=_t && j<_n+_t );
	 return _data[(i-_s) + (j-_t)*_m];
  }
  F& operator()(int i, int j)  {	 assert( i>=_s && i<_m+_s && j>=_t && j<_n+_t );
	 return _data[(i-_s) + (j-_t)*_m];
  }
  int m() const { return _m; }
  int n() const { return _n; }
  int s() const { return _s; }
  int t() const { return _t; }
  F* data() const { return _data; }
};

//INPUT template <class F> inline istream& operator>>(istream& is, OffMat<F>& vec)
//OUTPUT template <class F> inline ostream& operator<<( ostream& os, const OffMat<F>& mat)

template <class F> inline void setvalue(OffMat<F>& mat, F val)
{
  for(int i=mat.s(); i<mat.s()+mat.m(); i++)
	for(int j=mat.t(); j<mat.t()+mat.n(); j++)
	  mat(i,j) = val;
}

//template <class F> inline void clear(OffMat<F>& M){  memset(M.data(), 0, M.m()*M.n()*sizeof(F));}

typedef OffMat<bool>   BolOffMat;
typedef OffMat<char>   ChrOffMat;
typedef OffMat<int>    IntOffMat;
typedef OffMat<double> DblOffMat;
typedef OffMat<cpx>    CpxOffMat;


#endif

