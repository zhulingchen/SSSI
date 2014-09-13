/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#ifndef _OFFTNS_HPP_
#define _OFFTNS_HPP_

#include "commoninc.hpp"

//-------------------------------------------------
template <class F>
class OffTns {
public:
  int _m, _n, _p;
  int _s, _t, _u;
  F* _data;
public:
  OffTns(int m=0, int n=0, int p=0): _m(m), _n(n), _p(p), _s(-m/2), _t(-n/2), _u(-p/2) {
	 if(_m>0 && _n>0 && _p>0) {		_data = new F[_m*_n*_p]; assert( _data!=NULL );	 memset(_data, 0, _m*_n*_p*sizeof(F));	 } else		_data = NULL;
  }
  OffTns(const OffTns& C): _m(C._m), _n(C._n), _p(C._p), _s(C._s), _t(C._t), _u(C._u) {
	 if(_m>0 && _n>0 && _p>0) {		_data = new F[_m*_n*_p]; assert( _data!=NULL );	 memset(_data, 0, _m*_n*_p*sizeof(F));	 } else 		_data = NULL;
	 if(_m>0 && _n>0 && _p>0) {		memcpy( _data, C._data, _m*_n*_p*sizeof(F) );	 }
  }
  ~OffTns() {
	 if(_m>0 && _n>0 && _p>0) {		delete[] _data; _data = NULL; }
  }
  OffTns& operator=(const OffTns& C) {
	 if(_m>0 && _n>0 && _p>0) {		delete[] _data; _data = NULL; }
	 _m=C._m; _n=C._n; _p=C._p; _s=C._s; _t=C._t; _u=C._u;
	 if(_m>0 && _n>0 && _p>0) {		_data = new F[_m*_n*_p]; assert( _data!=NULL );	 memset(_data, 0, _m*_n*_p*sizeof(F));	 } else		_data = NULL;
	 if(_m>0 && _n>0 && _p>0) {		memcpy( _data, C._data, _m*_n*_p*sizeof(F) );	 }
	 return *this;
  }
  void resize(int m, int n, int p) {
	 if(_m!=m || _n!=n || _p!=p) {
		if(_m>0 && _n>0 && _p>0) {		delete[] _data; _data = NULL;		}
		_m=m; _n=n; _p=p; _s=-m/2; _t=-n/2; _u=-p/2;
		if(_m>0 && _n>0 && _p>0) {		_data = new F[_m*_n*_p]; assert( _data!=NULL );		memset(_data, 0, _m*_n*_p*sizeof(F));		} else		  _data = NULL;
	 }
  }
  const F& operator()(int i, int j, int k) const {	 assert( i>=_s && i<_m+_s && j>=_t && j<_n+_t && k>=_u && k<_p+_u);
	 return _data[(i-_s) + (j-_t)*_m + (k-_u)*_m*_n];
  }
  F& operator()(int i, int j, int k) {	 assert( i>=_s && i<_m+_s && j>=_t && j<_n+_t && k>=_u && k<_p+_u);
	 return _data[(i-_s) + (j-_t)*_m + (k-_u)*_m*_n];
  }
  int m() const { return _m; }
  int n() const { return _n; }
  int p() const { return _p; }
  int s() const { return _s; }
  int t() const { return _t; }
  int u() const { return _u; }
  F* data() const { return _data; }
};

template <class F> inline void setvalue(OffTns<F>& T, F val)
{
  for(int i=T.s(); i<T.m()+T.s(); i++)
	 for(int j=T.t(); j<T.n()+T.t(); j++)
		for(int k=T.u(); k<T.p()+T.u(); k++)
		  T(i,j,k) = val;
}
//template <class F> inline void clear(OffTns<F>& T){  memset(T.data(), 0, T.m()*T.n()*T.p()*sizeof(F));}

typedef OffTns<bool>   BolOffTns;
typedef OffTns<char>   ChrOffTns;
typedef OffTns<int>    IntOffTns;
typedef OffTns<double> DblOffTns;
typedef OffTns<cpx>    CpxOffTns;

//-------------------------------------------------
inline double energy(CpxOffTns& t)
{
  double val=0;
  cpx* data = t.data();
  for(int i=0; i<t.m()*t.n()*t.p(); i++) {	 val += norm(data[i]);  }
  return val;
}

#endif
