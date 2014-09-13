/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#ifndef _CPXNUMTNSBLKD_HPP_
#define _CPXNUMTNSBLKD_HPP_

#include "numtns.hpp"

//-----------------------------------------
//Complex Numerical Tensor Blocked (in all directions)
class CpxNumTnsBlkd
{
protected:  //PaData _padata;
  //input
  int _m, _n, _p; //length in each direction
  int _b; //size of the block in each direction
  IntNumTns _owners; //for each block, the owner processor
  IntNumTns _sizes;  //for each block, the size (in terms of the number of doubles)
  BolNumTns _exists; //for each box, whether exisits on thiis processor
  NumTns<CpxNumTns> _blocks;  //NumTns<CpxNumTns> _blocks; //a tensor of tensors
  //local  //int _e, _f, _g; //number of block in each direciton
public:
  CpxNumTnsBlkd() {;} //empty constructor
  CpxNumTnsBlkd(const CpxNumTnsBlkd& D); //copy
  ~CpxNumTnsBlkd() {;}
  CpxNumTnsBlkd& operator=(const CpxNumTnsBlkd& D); //assignment
  
  int setup(int m, int n, int p, int b, IntNumTns& owners); //setup data
  int expand(BolNumTns& newexists); //no communication, just expand with zero data
  int scatter(BolNumTns& newexists); //scatter from the owner to other procs
  int shift(IntNumTns& newowners); //shift owners
  int discard(); //discard non-owned blocks
  int combine(); //combine information from procs to the owner
  //access data
  int m() { return _m; }  int n() { return _n; }  int p() { return _p; }
  int b() { return _b; }
  int e() { return _m/_b; }  int f() { return _n/_b; }  int g() { return _p/_b; }
  int numblocks() { return e()*f()*g(); }
  
  IntNumTns& owners() { return _owners; }
  IntNumTns& sizes()  { return _sizes; }
  BolNumTns& exists() { return _exists; }
  
  CpxNumTns& block(int i, int j, int k) {
	 assert(_exists(i,j,k)==true);
	 CpxNumTns& tmp = _blocks(i,j,k);	 assert(tmp.m()==_b && tmp.n()==_b && tmp.p()==_b);
	 return _blocks(i,j,k);
  }
  double globalenergy();
  int check();
  
  //extra
  int mpirank() const { int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); return rank; }
  int mpisize() const { int size; MPI_Comm_size(MPI_COMM_WORLD, &size); return size; }
};

#endif
