/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#ifndef _NUMTNS_HPP_
#define _NUMTNS_HPP_

#include "commoninc.hpp"

template <class F>
class NumTns {
public:
  int _m, _n, _p;
  F* _data;
public:
  NumTns(int m=0, int n=0, int p=0): _m(m), _n(n), _p(p) {
	 if(_m>0 && _n>0 && _p>0) {		_data = new F[_m*_n*_p]; assert( _data!=NULL );	 } else		_data = NULL;
  }
  NumTns(const NumTns& C): _m(C._m), _n(C._n), _p(C._p) {
	 if(_m>0 && _n>0 && _p>0) {		_data = new F[_m*_n*_p]; assert( _data!=NULL );	 } else 		_data = NULL;
	 if(_m>0 && _n>0 && _p>0) {		for(int i=0; i<_m*_n*_p; i++) _data[i] = C._data[i]; }
  }
  ~NumTns() {
	 if(_m>0 && _n>0 && _p>0) {		delete[] _data; _data = NULL; }
  }
  NumTns& operator=(const NumTns& C) {
	 if(_m>0 && _n>0 && _p>0) {		delete[] _data; _data = NULL; }
	 _m = C._m; _n=C._n; _p=C._p;
	 if(_m>0 && _n>0 && _p>0) {		_data = new F[_m*_n*_p]; assert( _data!=NULL );	 } else		_data = NULL;
	 if(_m>0 && _n>0 && _p>0) {		for(int i=0; i<_m*_n*_p; i++) _data[i] = C._data[i]; }
	 return *this;
  }
  void resize(int m, int n, int p) {
	 if(_m!=m || _n!=n || _p!=p) {
		if(_m>0 && _n>0 && _p>0) {		delete[] _data; _data = NULL;		}
		_m = m; _n = n; _p=p;
		if(_m>0 && _n>0 && _p>0) {		_data = new F[_m*_n*_p]; assert( _data!=NULL );		} else		  _data = NULL;
	 }
  }
  const F& operator()(int i, int j, int k) const {	 assert( i>=0 && i<_m && j>=0 && j<_n && k>=0 && k<_p );
	 return _data[i + j*_m + k*_m*_n];
  }
  F& operator()(int i, int j, int k) {	 assert( i>=0 && i<_m && j>=0 && j<_n && k>=0 && k<_p );
	 return _data[i + j*_m + k*_m*_n];
  }
  int m() const { return _m; }
  int n() const { return _n; }
  int p() const { return _p; }
  F* data() const { return _data; }
};

template <class F> inline void setvalue(NumTns<F>& T, F val)
{
  for(int i=0; i<T.m(); i++)
	for(int j=0; j<T.n(); j++)
	  for(int k=0; k<T.p(); k++)
		T(i,j,k) = val;
}
//template <class F> inline void clear(NumTns<F>& T){  memset(T.data(), 0, T.m()*T.n()*T.p()*sizeof(F));}

//template <class F> inline istream& operator>>(istream& is, NumTns<F>& tns);
//template <class F> inline ostrema& operator<<(ostream& os, const NumTns<F>& tns);

typedef NumTns<bool>   BolNumTns;
typedef NumTns<char>   ChrNumTns;
typedef NumTns<int>    IntNumTns;
typedef NumTns<double> DblNumTns;
typedef NumTns<cpx>    CpxNumTns;

//-------------------------------------------------
inline double energy(CpxNumTns& t)
{
  double val=0;
  cpx* data = t.data();
  for(int i=0; i<t.m()*t.n()*t.p(); i++) {	 val += norm(data[i]);  }
  return val;
}


#endif
