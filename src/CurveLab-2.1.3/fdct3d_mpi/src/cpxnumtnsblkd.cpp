/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#include "cpxnumtnsblkd.hpp"

//---------------------------------------------
CpxNumTnsBlkd::CpxNumTnsBlkd(const CpxNumTnsBlkd& D):
  _m(D._m), _n(D._n), _p(D._p), _b(D._b), _owners(D._owners), _sizes(D._sizes), _exists(D._exists), _blocks(D._blocks)
{
}

CpxNumTnsBlkd& CpxNumTnsBlkd::operator=(const CpxNumTnsBlkd& D)
{
  _m = D._m;  _n = D._n;  _p = D._p;  _b = D._b;
  _owners = D._owners;
  _sizes = D._sizes;
  _exists = D._exists;
  _blocks = D._blocks;
  return *this;
}

//---------------------------------------------
double CpxNumTnsBlkd::globalenergy()
{
  double lclsum = 0;
  for(int i=0; i<e(); i++)
	 for(int j=0; j<f(); j++)
		for(int k=0; k<g(); k++)
		  if(_owners(i,j,k)==mpirank())
			 lclsum += energy( _blocks(i,j,k) );
  double glbsum = 0;
  iC( MPI_Reduce((void*)(&lclsum), (void*)(&glbsum), 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD) );
  return glbsum;
}

//---------------------------------------------
int CpxNumTnsBlkd::check()
{
  for(int i=0; i<e(); i++)	 for(int j=0; j<f(); j++)		for(int k=0; k<g(); k++)
	 if(_exists(i,j,k)==true) {
		CpxNumTns& tmp = _blocks(i,j,k);
		iA( tmp.m()==_b && tmp.n()==_b && tmp.p()==_b );
	 }
  return 0;
}

//---------------------------------------------
int CpxNumTnsBlkd::setup(int m, int n, int p, int b, IntNumTns& owners)
{
  _m = m;  _n = n;  _p = p;  _b = b;
  //owner
  iA( owners.m()==e() && owners.n()==f() && owners.p()==g() );
  _owners = owners;
  //sizes
  _sizes.resize(e(),f(),g());
  int C2D = 2;
  int blksz = _b*_b*_b * C2D;
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++)
	 _sizes(i,j,k) = blksz;
  //exists + blocks
  _exists.resize(e(),f(),g());
  _blocks.resize(e(),f(),g()); //all blocks are empty now
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++)
	 if(_owners(i,j,k)==mpirank()) {
		_exists(i,j,k) = true;
		_blocks(i,j,k).resize(_b,_b,_b);
	 } else {
		_exists(i,j,k) = false;
		_blocks(i,j,k).resize(0,0,0);
	 }
  return 0;
}

int CpxNumTnsBlkd::expand(BolNumTns& newexists)
{
  //LEXING: usually only called once
  //the final exists is the union of _exists and newexists
  iA( newexists.m()==e() && newexists.n()==f() && newexists.p()==g() );
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++)
	 if(_exists(i,j,k)==false && newexists(i,j,k)==true) {
		_exists(i,j,k) = true;
		_blocks(i,j,k).resize(_b,_b,_b);
	 }
  return 0;
}

int CpxNumTnsBlkd::scatter(BolNumTns& newexists)
{
  //LEXING: usually only called once
  //1. the global vector
  vector<int> glblszs(mpisize(), 0);
  int glbnum = 0;
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++) {
	 int pi = _owners(i,j,k);
	 glblszs[pi] += _sizes(i,j,k);
	 glbnum += _sizes(i,j,k);
  }
  
  vector<int> glbaccs(mpisize(), 0);
  int tmp = 0;
  for(int pi=0; pi<mpisize(); pi++) {
	 glbaccs[pi] = tmp;
	 tmp += glblszs[pi];
  }
  IntNumTns glbstts(e(),f(),g());
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++) {
	 int pi = _owners(i,j,k);
	 glbstts(i,j,k) = glbaccs[pi];
	 glbaccs[pi] += _sizes(i,j,k);
  }
  
  int lclsum = 0;
  vector<int> lclbid;
  vector<int> glbbid;
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++) {
	 if(newexists(i,j,k)==true && _exists(i,j,k)==false) {
		lclbid.push_back(lclsum);
		glbbid.push_back(glbstts(i,j,k));
		lclsum += _sizes(i,j,k);
	 }
  }
  int C2D = 2;
  int blksz = _b*_b*_b * C2D;
  IS lclis;  iC( ISCreateBlock(PETSC_COMM_SELF,  blksz, lclbid.size(), &(lclbid[0]), &lclis) );
  IS glbis;  iC( ISCreateBlock(PETSC_COMM_WORLD, blksz, glbbid.size(), &(glbbid[0]), &glbis) );
  lclbid.clear();
  glbbid.clear();

  /*
  vector<PetscInt> l2gmap;
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++) {
	 if(newexists(i,j,k)==true && _exists(i,j,k)==false) {
		lclsum += _sizes(i,j,k);
		for(int g=0; g<_sizes(i,j,k); g++)
		  l2gmap.push_back( glbstts(i,j,k) + g );
	 }
  }
  iA(l2gmap.size()==lclsum);
  IS lclis;  iC( ISCreateStride(PETSC_COMM_SELF, l2gmap.size(), 0, 1, &lclis) );
  IS glbis;  iC( ISCreateGeneral(PETSC_COMM_WORLD, l2gmap.size(), &(l2gmap[0]), &glbis) );
  l2gmap.clear(); //SAVE SPACE
  */
  
  //2. allocate a global vector, and copy data
  Vec glbvec;  iC( VecCreateMPI(PETSC_COMM_WORLD, glblszs[mpirank()], PETSC_DETERMINE, &glbvec) );
  double* glbarr;  iC( VecGetArray(glbvec, &glbarr) );
  double* glbptr = glbarr;
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++) {
	 int pi = _owners(i,j,k);
	 if(pi==mpirank()) {
		double* tmpptr = (double*)(_blocks(i,j,k).data());
		for(int g=0; g<_sizes(i,j,k); g++) {
		  *glbptr = tmpptr[g];		  glbptr++;
		}
	 }
  }
  iC( VecRestoreArray(glbvec, &glbarr) );
  
  Vec lclvec;  iC( VecCreateSeq(PETSC_COMM_SELF, lclsum, &lclvec) );  //cerr<<lclsum<<endl;
  
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
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++) {
	 if(newexists(i,j,k)==true && _exists(i,j,k)==false) {
		_blocks(i,j,k).resize(_b,_b,_b);
		double* tmpptr = (double*)(_blocks(i,j,k).data());
		for(int g=0; g<_sizes(i,j,k); g++) {
		  tmpptr[g] = *lclptr;		  lclptr++;
		}
		_exists(i,j,k) = true; //VERY IMPORTANT
	 }
  }
  iC( VecRestoreArray(lclvec, &lclarr) );
  iC( VecDestroy(lclvec) );
  return 0;
}

int CpxNumTnsBlkd::shift(IntNumTns& newowners)
{
  assert( newowners.m()==e() && newowners.n()==f() && newowners.p()==g() );
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++) {
	if(newowners(i,j,k)==mpirank()) {
	  iA(_exists(i,j,k)==true);
	}
  }
  _owners = newowners;
  return 0;
}

int CpxNumTnsBlkd::discard()
{
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++) {
	if(_owners(i,j,k)!=mpirank()) {
	  _blocks(i,j,k).resize(0,0,0);
	  _exists(i,j,k) = false;
	}
  }
  return 0;
}

int CpxNumTnsBlkd::combine()
{
  //LEXING: usually only called once
  //1. the global vector
  vector<int> glblszs(mpisize(), 0);
  int glbnum = 0;
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++) {
	 int pi = _owners(i,j,k);
	 glblszs[pi] += _sizes(i,j,k);
	 glbnum += _sizes(i,j,k);
  }  //iC( MPI_Barrier(MPI_COMM_WORLD) );  iC( PetscPrintf(MPI_COMM_WORLD, "combine glbnum %d\n", glbnum) );
  
  vector<int> glbaccs(mpisize(), 0);
  int tmp = 0;
  for(int pi=0; pi<mpisize(); pi++) {
	 glbaccs[pi] = tmp;
	 tmp += glblszs[pi];
  }
  IntNumTns glbstts(e(),f(),g());
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++) {
	 int pi = _owners(i,j,k);
	 glbstts(i,j,k) = glbaccs[pi];
	 glbaccs[pi] += _sizes(i,j,k);
  }
  
  int lclsum = 0;
  vector<int> lclbid;
  vector<int> glbbid;
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++) {
	 if(_exists(i,j,k)==true && _owners(i,j,k)!=mpirank()) {
		lclbid.push_back(lclsum);
		glbbid.push_back(glbstts(i,j,k));
		lclsum += _sizes(i,j,k);
	 }
  }
  int C2D = 2;
  int blksz = _b*_b*_b * C2D;
  IS lclis;  iC( ISCreateBlock(PETSC_COMM_SELF,  blksz, lclbid.size(), &(lclbid[0]), &lclis) );
  IS glbis;  iC( ISCreateBlock(PETSC_COMM_WORLD, blksz, glbbid.size(), &(glbbid[0]), &glbis) );
  lclbid.clear();
  glbbid.clear();

  /*
  int lclsum = 0;
  vector<PetscInt> l2gmap;
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++) {
	 if(_exists(i,j,k)==true && _owners(i,j,k)!=mpirank()) {
		lclsum += _sizes(i,j,k);
		for(int g=0; g<_sizes(i,j,k); g++)
		  l2gmap.push_back( glbstts(i,j,k) + g );
	 }
  }
  iA(l2gmap.size()==lclsum);
  IS lclis;  iC( ISCreateStride(PETSC_COMM_SELF, l2gmap.size(), 0, 1, &lclis) );
  IS glbis;  iC( ISCreateGeneral(PETSC_COMM_WORLD, l2gmap.size(), &(l2gmap[0]), &glbis) );
  l2gmap.clear(); //SAVE SPACE
  */  

  //2. allocate a global vector and a local vector, put data in local
  Vec glbvec;  iC( VecCreateMPI(PETSC_COMM_WORLD, glblszs[mpirank()], PETSC_DETERMINE, &glbvec) );

  Vec lclvec;  iC( VecCreateSeq(PETSC_COMM_SELF, lclsum, &lclvec) );  //cerr<<lclsum<<endl;
  double* lclarr;  iC( VecGetArray(lclvec, &lclarr) );
  double* lclptr = lclarr;
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++) {
	 if(_exists(i,j,k)==true && _owners(i,j,k)!=mpirank()) {
		double* tmpptr = (double*)(_blocks(i,j,k).data());
		for(int g=0; g<_sizes(i,j,k); g++) {
		  *lclptr = tmpptr[g];		  lclptr++;
		}
	 }
  }
  iC( VecRestoreArray(lclvec, &lclarr) );
  
  //3. vec scatter
  VecScatter sc;  iC( VecScatterCreate(glbvec, glbis, lclvec, lclis, &sc) );
  iC( ISDestroy(lclis) );  iC( ISDestroy(glbis) ); //SAVE SPACE
  
  iC( VecScatterBegin(lclvec, glbvec, ADD_VALUES, SCATTER_REVERSE, sc) );
  iC( VecScatterEnd(  lclvec, glbvec, ADD_VALUES, SCATTER_REVERSE, sc) );
  
  iC( VecScatterDestroy(sc) ); //SAVE SPACE
  iC( VecDestroy(lclvec) );
  
  //4. store
  double* glbarr;  iC( VecGetArray(glbvec, &glbarr) );
  double* glbptr = glbarr;
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++) {
	 int pi = _owners(i,j,k);
	 if(pi==mpirank()) {
		double* tmpptr = (double*)(_blocks(i,j,k).data());
		for(int g=0; g<_sizes(i,j,k); g++) {
		  tmpptr[g] += *glbptr;		  glbptr++; //LEXING: += is very important
		}
	 }
  }
  iC( VecRestoreArray(glbvec, &glbarr) );
  iC( VecDestroy(glbvec) );
  
  //IMPORTANT
  for(int k=0; k<g(); k++)	 for(int j=0; j<f(); j++)		for(int i=0; i<e(); i++) {
	if(_owners(i,j,k)!=mpirank()) {
	  _blocks(i,j,k).resize(0,0,0);
	  _exists(i,j,k) = false;
	}
  }
  
  return 0;
}
