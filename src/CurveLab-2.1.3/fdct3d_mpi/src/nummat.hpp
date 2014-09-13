/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#ifndef _NUMMAT_HPP_
#define _NUMMAT_HPP_

#include "commoninc.hpp"

//-------------------------------------------------
template <class F>
class NumMat {
public:
  int _m, _n;
  F* _data;
public:
  NumMat(int m=0, int n=0): _m(m), _n(n) {
	 if(_m>0 && _n>0) {		_data = new F[_m*_n]; assert( _data!=NULL );	 memset(_data, 0, _m*_n*sizeof(F));	 } else		_data = NULL;
  }
  NumMat(const NumMat& C): _m(C._m), _n(C._n) {
	 if(_m>0 && _n>0) {		_data = new F[_m*_n]; assert( _data!=NULL );	 memset(_data, 0, _m*_n*sizeof(F));	 } else 		_data = NULL;
	 if(_m>0 && _n>0) {		memcpy( _data, C._data, _m*_n*sizeof(F) );	 }
  }
  ~NumMat() {
	 if(_m>0 && _n>0) {		delete[] _data; _data = NULL; }
  }
  NumMat& operator=(const NumMat& C) {
	 if(_m>0 && _n>0) {		delete[] _data; _data = NULL; }
	 _m = C._m; _n=C._n;
	 if(_m>0 && _n>0) {		_data = new F[_m*_n]; assert( _data!=NULL );	 memset(_data, 0, _m*_n*sizeof(F));	 } else		_data = NULL;
	 if(_m>0 && _n>0) {		memcpy( _data, C._data, _m*_n*sizeof(F) );	 }
	 return *this;
  }
  void resize(int m, int n)  {
	 if(_m!=m || _n!=n) {
		if(_m>0 && _n>0) {		delete[] _data; _data = NULL;		}
		_m = m; _n = n;
		if(_m>0 && _n>0) {		_data = new F[_m*_n]; assert( _data!=NULL );		memset(_data, 0, _m*_n*sizeof(F));		} else		  _data = NULL;
	 }
  }
  const F& operator()(int i, int j) const  { 	 assert( i>=0 && i<_m && j>=0 && j<_n );
	 return _data[i + j*_m];
  }
  F& operator()(int i, int j)  { 	 assert( i>=0 && i<_m && j>=0 && j<_n );
	 return _data[i + j*_m];
  }
  int m() const { return _m; }
  int n() const { return _n; }
  F* data() const { return _data; }
};

//INPUT template <class F> inline istream& operator>>(istream& is, NumMat<F>& mat)
//OUTPUT template <class F> inline ostream& operator<<(ostream& os, const NumMat<F>& mat)

template <class F> inline void setvalue(NumMat<F>& M, F val)
{
  for(int i=0; i<M.m(); i++)
	 for(int j=0; j<M.n(); j++)
		M(i,j) = val;
}
//template <class F> inline void clear(NumMat<F>& M){  memset(M.data(), 0, M.m()*M.n()*sizeof(F));}

typedef NumMat<bool>   BolNumMat;
typedef NumMat<char>   ChrNumMat;
typedef NumMat<int>    IntNumMat;
typedef NumMat<double> DblNumMat;
typedef NumMat<cpx>    CpxNumMat;

#endif

