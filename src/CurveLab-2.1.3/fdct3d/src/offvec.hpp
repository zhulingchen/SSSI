/*	FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#ifndef _OFFVEC_HPP_
#define _OFFVEC_HPP_

#include "commoninc.hpp"

//-------------------------------------------------
template <class F>
class OffVec {
public:
  int _m;
  int _s; //starting position, using negative
  F* _data;
public:
  OffVec(int m=0): _m(m), _s(-m/2) {
	 if(_m>0) {		_data = new F[_m]; assert( _data!=NULL );	 memset(_data, 0, _m*sizeof(F));	 } else		_data=NULL;
  }
  OffVec(const OffVec& C): _m(C._m), _s(C._s) {
	 if(_m>0) {		_data = new F[_m]; assert( _data!=NULL );	 memset(_data, 0, _m*sizeof(F));	 } else 		_data=NULL;
	 if(_m>0) {		memcpy( _data, C._data, _m*sizeof(F) );	 }
  }
  ~OffVec() {
	 if(_m>0) {		delete[] _data; _data = NULL;	 }
  }
  OffVec& operator=(const OffVec& C) {
	 if(_m>0) {		delete[] _data; _data = NULL;	 }
	 _m = C._m;	 _s = C._s;
	 if(_m>0) {		_data = new F[_m]; assert( _data!=NULL );	 memset(_data, 0, _m*sizeof(F));	 } else		_data = NULL;
	 if(_m>0) {		memcpy( _data, C._data, _m*sizeof(F) );	 }
	 return *this;
  }
  void resize(int m)  {
	 if(_m!=m) {
		if(_m>0) {		  delete[] _data; _data = NULL;		}
		_m = m; _s = -m/2;
		if(_m>0) {		  _data = new F[_m]; assert( _data!=NULL );		memset(_data, 0, _m*sizeof(F));		}		else		  _data = NULL;
	 }
  }
  const F& operator()(int i) const  { 	 //assert( i>=_s && i<_m+_s);
	 return _data[i-_s];
  }
  F& operator()(int i)  { 	 //assert( i>=_s && i<_m+_s);
	 return _data[i-_s];
  }
  int m() const { return _m; }
  int s() const { return _s; }
  F* data() const { return _data; }
};

/*
//INPUT
template <class F> inline istream& operator>>(istream& is, OffVec<F>& vec)
{
  int m; is>>m;
  vec.resize(m);
  for(int i=0; i<vec.m(); i++)
	 is>>vec(i);
  return is;
}
*/

//OUTPUT
template <class F> inline ostream& operator<<(ostream& os, const OffVec<F>& vec)
{
  os<<vec.m()<<" "<<vec.s()<<endl;
  os.setf(ios_base::scientific, ios_base::floatfield);
  for(int i=vec.s(); i<vec.s()+vec.m(); i++)
	 os<<" "<<vec(i)<<endl;
  return os;
}

//SET VALUE
template <class F> inline void setvalue(OffVec<F>& vec, F val)
{
  for(int i=vec.s(); i<vec.s()+vec.m(); i++)
	 vec(i) = val;
}
//CLEAR
template <class F> inline void clear(OffVec<F>& vec)
{
  memset(vec.data(), 0, vec.m()*sizeof(F));
}

typedef OffVec<int>    IntOffVec;
typedef OffVec<double> DblOffVec;
typedef OffVec<cpx>    CpxOffVec;

#endif

