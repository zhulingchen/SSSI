/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#include "cpxcrvletprtd.hpp"

CpxCrvletPrtd::CpxCrvletPrtd(const CpxCrvletPrtd& D):
  _nx(D._nx), _ny(D._ny), _nz(D._nz), _owners(D._owners), _sizes(D._sizes), _exists(D._exists), _blocks(D._blocks)
{}

CpxCrvletPrtd& CpxCrvletPrtd::operator=(const CpxCrvletPrtd& D)
{
  _nx = D._nx;  _ny = D._ny;  _nz = D._nz;
  _owners = D._owners;
  _sizes = D._sizes;
  _exists = D._exists;
  _blocks = D._blocks;
  return *this;
}

//---------------------------------------------
double CpxCrvletPrtd::globalenergy()
{
  double lclsum = 0;
  vector< vector<int> >& c = _nx;
  for(int s=0; s<c.size(); s++)
	 for(int w=0; w<c[s].size(); w++)
		if(_owners[s][w]==mpirank())
		  lclsum += energy(_blocks[s][w]);
  double glbsum = 0;
  iC( MPI_Reduce((void*)(&lclsum), (void*)(&glbsum), 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD) );
  return glbsum;
}
//---------------------------------------------
int CpxCrvletPrtd::check()
{
  vector< vector<int> >& c = _nx;
  for(int s=0; s<c.size(); s++)
	 for(int w=0; w<c[s].size(); w++)
		if(_exists[s][w]==true) {
		  CpxNumTns& tmp = _blocks[s][w];
		  iA( tmp.m()==_nx[s][w] && tmp.n()==_ny[s][w] && tmp.p()==_nz[s][w] );
		}
  return 0;
}

//-------------------------------------------------
int CpxCrvletPrtd::setup(vector< vector<int> > nx, vector< vector<int> > ny, vector< vector<int> > nz, vector< vector<int> >& owners)
{
  _nx = nx;  _ny = ny;  _nz = nz;
  
  //owners
  //check owners's size is the same as nx
  _owners = owners;
  //sizes
  int C2D = 2;
  vector< vector<int> >& c = _nx;
  _sizes.resize(nx.size());
  for(int s=0; s<c.size(); s++) {	 _sizes[s].resize( c[s].size() );
	 for(int w=0; w<c[s].size(); w++)
		_sizes[s][w] = _nx[s][w]*_ny[s][w]*_nz[s][w] * C2D;
  }
  //exists + blocks
  _exists.resize(c.size());
  _blocks.resize(c.size());
  for(int s=0; s<c.size(); s++) {	 _exists[s].resize( c[s].size() );	 _blocks[s].resize( c[s].size() );
	 for(int w=0; w<c[s].size(); w++)
		if(_owners[s][w]==mpirank()) {
		  _exists[s][w] = true;
		  _blocks[s][w].resize(_nx[s][w], _ny[s][w], _nz[s][w]);
		} else {
		  _exists[s][w] = false;
		  _blocks[s][w].resize(0,0,0);
		}
  }
  return 0;
}

int CpxCrvletPrtd::expand(vector< vector<bool> >& newexists)
{
  //todo: check size of newexists is right
  //the final exists is the union of _exists and newexisits
  vector< vector<int> >& c = _nx;
  for(int s=0; s<c.size(); s++) {
	 for(int w=0; w<c[s].size(); w++)
		if(_exists[s][w]==false && newexists[s][w]==true) {
		  _exists[s][w] = true;
		  _blocks[s][w].resize(_nx[s][w], _ny[s][w], _nz[s][w]);
		}
  }
  return 0;
}

int CpxCrvletPrtd::scatter(vector< vector<bool> >& newexists)
{
  //LEXING: usually only called once
  vector< vector<int> >& c = _nx;
  //1. the global vector
  vector<int> glblszs(mpisize(), 0);
  int glbnum = 0;
  for(int s=0; s<c.size(); s++)	 for(int w=0; w<c[s].size(); w++) {
	 int pi = _owners[s][w];
	 glblszs[pi] += _sizes[s][w];
	 glbnum += _sizes[s][w];
  }
  vector<int> glbaccs(mpisize(), 0);
  int tmp = 0;
  for(int pi=0; pi<mpisize(); pi++) {
	 glbaccs[pi] = tmp;
	 tmp += glblszs[pi];
  }
  vector< vector<int> > glbstts(c); //not cleared, but okay
  for(int s=0; s<c.size(); s++)	 for(int w=0; w<c[s].size(); w++) {
	 int pi = _owners[s][w];
	 glbstts[s][w] = glbaccs[pi];
	 glbaccs[pi] += _sizes[s][w];
  }
  
  int lclsum = 0;
  vector<int> l2gmap;
  for(int s=0; s<c.size(); s++)	 for(int w=0; w<c[s].size(); w++) {
	 if(newexists[s][w]==true && _exists[s][w]==false) {
		lclsum += _sizes[s][w];
		for(int g=0; g<_sizes[s][w]; g++)
		  l2gmap.push_back( glbstts[s][w] + g );
	 }
  }
  iA(l2gmap.size()==lclsum);
  
  IS lclis;  iC( ISCreateStride(PETSC_COMM_SELF, l2gmap.size(), 0, 1, &lclis) );
  IS glbis;  iC( ISCreateGeneral(PETSC_COMM_WORLD, l2gmap.size(), &(l2gmap[0]), &glbis) );
  l2gmap.clear(); //SAVE SPACE
  
  //2. allocate a global vector, and copy data
  Vec glbvec;  iC( VecCreateMPI(PETSC_COMM_WORLD, glblszs[mpirank()], PETSC_DETERMINE, &glbvec) );
  double* glbarr;  iC( VecGetArray(glbvec, &glbarr) );
  double* glbptr = glbarr;
  for(int s=0; s<c.size(); s++)	 for(int w=0; w<c[s].size(); w++) {
	int pi = _owners[s][w];
	if(pi==mpirank()) {
	  double* tmpptr = (double*)(_blocks[s][w].data());
	  for(int g=0; g<_sizes[s][w]; g++) {
		*glbptr = tmpptr[g];		  glbptr++;
	  }
	}
  }
  iC( VecRestoreArray(glbvec, &glbarr) );
  
  Vec lclvec;  iC( VecCreateSeq(PETSC_COMM_SELF, lclsum, &lclvec) );
  
  //3. vec scatter
  VecScatter sc;  iC( VecScatterCreate(glbvec, glbis, lclvec, lclis, &sc) );
  iC( ISDestroy(lclis) );  iC( ISDestroy(glbis) ); //SAVE SPACE
  
  iC( VecScatterBegin(glbvec, lclvec, INSERT_VALUES, SCATTER_FORWARD, sc) );
  iC( VecScatterEnd(  glbvec, lclvec, INSERT_VALUES, SCATTER_FORWARD, sc) );
  
  iC( VecScatterDestroy(sc) ); //SAVE SPACE
  iC( VecDestroy(glbvec) );
  
  //4. store
  double* lclarr;  iC( VecGetArray(lclvec, &lclarr) );
  double* lclptr = lclarr;
  for(int s=0; s<c.size(); s++)	 for(int w=0; w<c[s].size(); w++) {
	 if(newexists[s][w]==true && _exists[s][w]==false) {
		_blocks[s][w].resize(_nx[s][w], _ny[s][w], _nz[s][w]);
		double* tmpptr = (double*)(_blocks[s][w].data());
		for(int g=0; g<_sizes[s][w]; g++) {
		  tmpptr[g] = *lclptr;		  lclptr++;
		}
		_exists[s][w] = true; //VERY IMPORTANT
	 }
  }
  iC( VecRestoreArray(lclvec, &lclarr) );
  iC( VecDestroy(lclvec) );
  return 0;
}

int CpxCrvletPrtd::shift(vector< vector<int> >& newowners)
{
  vector< vector<int> >& c = _nx;
  for(int s=0; s<c.size(); s++)	 for(int w=0; w<c[s].size(); w++) {
	 if(newowners[s][w]==mpirank()) {
		iA(_exists[s][w]==true);
	 }
  }
  _owners = newowners;
  return 0;
}

int CpxCrvletPrtd::discard()
{
  vector< vector<int> >& c = _nx;
  for(int s=0; s<c.size(); s++)	 for(int w=0; w<c[s].size(); w++) {
	 if(_owners[s][w]!=mpirank()) {
		_blocks[s][w].resize(0,0,0);
		_exists[s][w] = false;
	 }
  }
  return 0;
}

int CpxCrvletPrtd::combine()
{
  //LEXING: usually only called once
  vector< vector<int> >& c = _nx;
  //1. the global vector
  vector<int> glblszs(mpisize(), 0);
  int glbnum = 0;
  for(int s=0; s<c.size(); s++)	 for(int w=0; w<c[s].size(); w++) {
	 int pi = _owners[s][w];
	 glblszs[pi] += _sizes[s][w];
	 glbnum += _sizes[s][w];
  }
  vector<int> glbaccs(mpisize(), 0);
  int tmp = 0;
  for(int pi=0; pi<mpisize(); pi++) {
	 glbaccs[pi] = tmp;
	 tmp += glblszs[pi];
  }
  vector< vector<int> > glbstts(c);
  for(int s=0; s<c.size(); s++)	 for(int w=0; w<c[s].size(); w++) {
	 int pi = _owners[s][w];
	 glbstts[s][w] = glbaccs[pi];
	 glbaccs[pi] += _sizes[s][w];
  }
  
  int lclsum = 0;
  vector<int> l2gmap;
  for(int s=0; s<c.size(); s++)	 for(int w=0; w<c[s].size(); w++) {
	 if(_exists[s][w]==true && _owners[s][w]!=mpirank()) {
		lclsum += _sizes[s][w];
		for(int g=0; g<_sizes[s][w]; g++)
		  l2gmap.push_back( glbstts[s][w] + g );
	 }
  }
  iA(l2gmap.size()==lclsum);
  
  IS lclis;  iC( ISCreateStride(PETSC_COMM_SELF, l2gmap.size(), 0, 1, &lclis) );
  IS glbis;  iC( ISCreateGeneral(PETSC_COMM_WORLD, l2gmap.size(), &(l2gmap[0]), &glbis) );
  l2gmap.clear(); //SAVE SPACE
  
  //2. allocate a global vector and a local vector, put data in local
  Vec glbvec;  iC( VecCreateMPI(PETSC_COMM_WORLD, glblszs[mpirank()], PETSC_DETERMINE, &glbvec) );

  Vec lclvec;  iC( VecCreateSeq(PETSC_COMM_SELF, lclsum, &lclvec) );
  double* lclarr;  iC( VecGetArray(lclvec, &lclarr) );
  double* lclptr = lclarr;
  for(int s=0; s<c.size(); s++)	 for(int w=0; w<c[s].size(); w++) {
	 if(_exists[s][w]==true && _owners[s][w]!=mpirank()) {
		double* tmpptr = (double*)(_blocks[s][w].data());
		for(int g=0; g<_sizes[s][w]; g++) {
		  *lclptr = tmpptr[g];		  lclptr++;
		}
	 }
  }
  iC( VecRestoreArray(lclvec, &lclarr) );
  
  //3. vec scatter
  VecScatter sc;  iC( VecScatterCreate(glbvec, glbis, lclvec, lclis, &sc) );
  iC( ISDestroy(lclis) );  iC( ISDestroy(glbis) ); //SAVE SPACE
  
  iC( VecScatterBegin(glbvec, lclvec, ADD_VALUES, SCATTER_REVERSE, sc) );
  iC( VecScatterEnd(  glbvec, lclvec, ADD_VALUES, SCATTER_REVERSE, sc) );
  
  iC( VecScatterDestroy(sc) ); //SAVE SPACE
  iC( VecDestroy(lclvec) );
  
  //4. store
  double* glbarr;  iC( VecGetArray(glbvec, &glbarr) );
  double* glbptr = glbarr;
  for(int s=0; s<c.size(); s++)	 for(int w=0; w<c[s].size(); w++) {
	 int pi = _owners[s][w];
	 if(pi==mpirank()) {
		double* tmpptr = (double*)(_blocks[s][w].data());
		for(int g=0; g<_sizes[s][w]; g++) {
		  tmpptr[g] += *glbptr;		  glbptr++; //LEXING: += is very important
		}
	 }
  }
  iC( VecRestoreArray(glbvec, &glbarr) );
  iC( VecDestroy(glbvec) );
  
  //IMPORTANT
  for(int s=0; s<c.size(); s++)	 for(int w=0; w<c[s].size(); w++) {
	 if(_owners[s][w]!=mpirank()) {
		_blocks[s][w].resize(0,0,0);
		_exists[s][w] = false;
	 }
  }
  return 0;
}
