/*
   Copyright (C) 2004 Caltech
   Written by Lexing Ying
*/

#ifndef _NUMVEC_HPP_
#define _NUMVEC_HPP_

#include "fdct_wrapping_inc.hpp"

using std::ostream;
using std::istream;
using std::endl;

FDCT_WRAPPING_NS_BEGIN_NAMESPACE

//-------------------------------------------------
template <class F>
class NumVec {
public:
  int _m;
  F* _data;
public:
  NumVec(int m=0): _m(m) {
	 if(_m>0) {		_data = new F[_m]; assert( _data!=NULL );	 memset(_data, 0, _m*sizeof(F));	 } else		_data=NULL;
  }
  NumVec(const NumVec& C): _m(C._m) {
	 if(_m>0) {		_data = new F[_m]; assert( _data!=NULL );	 memset(_data, 0, _m*sizeof(F));	 } else 		_data=NULL;
	 if(_m>0) {		memcpy( _data, C._data, _m*sizeof(F) );	 }
  }
  ~NumVec() {
	 if(_m>0) {		delete[] _data; _data = NULL;	 }
  }
  NumVec& operator=(const NumVec& C) {
	 if(_m>0) {		delete[] _data; _data = NULL;	 }
	 _m = C._m;
	 if(_m>0) {		_data = new F[_m]; assert( _data!=NULL );	 memset(_data, 0, _m*sizeof(F));	 } else		_data = NULL;
	 if(_m>0) {		memcpy( _data, C._data, _m*sizeof(F) );	 }
	 return *this;
  }
  void resize(int m)  {
	 if(_m!=m) {
		if(_m>0) {		  delete[] _data; _data = NULL;		}
		_m = m;
		if(_m>0) {		  _data = new F[_m]; assert( _data!=NULL );		memset(_data, 0, _m*sizeof(F));		}		else		  _data = NULL;
	 }
  }
  const F& operator()(int i) const  {
 	 assert( i>=0 && i<_m);
	 return _data[i];
  }
  F& operator()(int i)  {
 	 assert( i>=0 && i<_m);
	 return _data[i];
  }
  int m() const { return _m; }
  F* data() const { return _data; }
};

//INPUT
template <class F> inline istream& operator>>(istream& is, NumVec<F>& vec)
{
  int m; is>>m;
  vec.resize(m);
  for(int i=0; i<vec.m(); i++)
	 is>>vec(i);
  return is;
}

//OUTPUT
template <class F> inline ostream& operator<<(ostream& os, const NumVec<F>& vec)
{
  os<<vec.m()<<endl;
  for(int i=0; i<vec.m(); i++)
	 os<<" "<<vec(i)<<endl;
  return os;
}

//SET VALUE
template <class F> inline void setvalue(NumVec<F>& vec, F val)
{
  for(int i=0; i<vec.m(); i++)
	vec(i) = val;
}
//CLEAR
template <class F> inline void clear(NumVec<F>& vec)
{
  memset(vec.data(), 0, vec.m()*sizeof(F));
}

typedef NumVec<int>    IntNumVec;
typedef NumVec<double> DblNumVec;
typedef NumVec<cpx>    CpxNumVec;

FDCT_WRAPPING_NS_END_NAMESPACE

#endif

