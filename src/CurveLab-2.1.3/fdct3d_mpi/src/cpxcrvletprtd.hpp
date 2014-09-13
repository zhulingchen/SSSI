/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#ifndef _CPXCRVLETPRTD_HPP_
#define _CPXCRVLETPRTD_HPP_

#include "numtns.hpp"

//-----------------------------------------
//Complex Curvelet Partitioned
class CpxCrvletPrtd
{
protected:
  vector< vector<int> > _nx, _ny, _nz; //size of
  vector< vector<int> > _owners;
  vector< vector<int> > _sizes;
  vector< vector<bool> > _exists;
  vector< vector<CpxNumTns> > _blocks;
  
public:
  CpxCrvletPrtd() {;}
  CpxCrvletPrtd(const CpxCrvletPrtd& D);
  ~CpxCrvletPrtd() {;}
  CpxCrvletPrtd& operator=(const CpxCrvletPrtd& D);
  int setup(vector< vector<int> > nx, vector< vector<int> > ny, vector< vector<int> > nz, vector< vector<int> >& owners);
  int expand(vector< vector<bool> >& newexists);
  int scatter(vector< vector<bool> >& newexists);
  int shift(vector< vector<int> >& newowners);
  int discard();
  int combine();
  //access
  vector< vector<int> >& nx() { return _nx; }
  vector< vector<int> >& ny() { return _ny; }
  vector< vector<int> >& nz() { return _nz; }
  
  vector< vector<int> >& owners() { return _owners; }
  vector< vector<int> >& sizes() { return _sizes; }
  vector< vector<bool> >& exists() { return _exists; }
  
  CpxNumTns& block(int s, int w) {	 assert(_exists[s][w]==true);	 return _blocks[s][w];  }
  double globalenergy();
  int check();
  
  //extra
  int mpirank() const { int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); return rank; }
  int mpisize() const { int size; MPI_Comm_size(MPI_COMM_WORLD, &size); return size; }
};

#endif
