/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#include "fdct3d.hpp"
#include "fdct3dinline.hpp"

int main(int argc, char** argv)
{
  PetscInitialize(&argc,&argv,"options",NULL);  //PetscTruth flg = PETSC_FALSE;
  srand48(0);
  
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  iC( PetscPrintf(MPI_COMM_WORLD, "mpisize %d\n", mpisize) );  
  
  PetscTruth flg = PETSC_FALSE;
  int m;  iC( PetscOptionsGetInt("", "-m", &m, &flg) ); iA(flg==PETSC_TRUE);
  int n;  iC( PetscOptionsGetInt("", "-n", &n, &flg) ); iA(flg==PETSC_TRUE);
  int p;  iC( PetscOptionsGetInt("", "-p", &p, &flg) ); iA(flg==PETSC_TRUE);
  int b;  iC( PetscOptionsGetInt("", "-b", &b, &flg) ); iA(flg==PETSC_TRUE);
  int nbscales;  iC( PetscOptionsGetInt("", "-nbscales", &nbscales, &flg) ); iA(flg==PETSC_TRUE);
  int nbdstz_coarse;  iC( PetscOptionsGetInt("", "-nbdstz_coarse", &nbdstz_coarse, &flg) ); iA(flg==PETSC_TRUE);
  
  CpxNumTnsBlkd X;
  BolNumTns newtnsexists(m/b,n/b,p/b);
  IntNumTns newtnsowners(m/b,n/b,p/b);
  iC( fdct3d_partition_cpxnumtnsblkd_z(m,n,p,b, newtnsexists, newtnsowners) );
  X.setup(m,n,p,b, newtnsowners);
  
  //1. generate data
  int e = X.e();  int f = X.f();  int g = X.g();
  for(int i=0; i<e; i++)	 for(int j=0; j<f; j++)		for(int k=0; k<g; k++) {
	 if(X.owners()(i,j,k)==mpirank) {
		CpxNumTns& Xblk = X.block(i,j,k);
		for(int ioff=0; ioff<b; ioff++)		  for(int joff=0; joff<b; joff++)			 for(int koff=0; koff<b; koff++) {
		  Xblk(ioff,joff,koff) = cpx(drand48(), 0); //cpx(drand48(), drand48());
		}
	 }
  }
  double xene = X.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "energy %e\n", xene) );  //  if(mpirank==0)	 cerr<<"X energy  "<<xene<<endl;
  
  time_t tm0, tm1;  tm0 = time(NULL);
  //2. forward
  CpxCrvletPrtd C;  CpxNumTnsBlkd W;
  iC( fdct3d_forward(m,n,p, nbscales,nbdstz_coarse, X, C,W) );
  tm1 = time(NULL);  iC( PetscPrintf(MPI_COMM_WORLD, "FORWARD %e\n", difftime(tm1,tm0)) );  tm0 = tm1;
  double cene = C.globalenergy();
  double wene = W.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "CW energy %e %e %e\n", cene, wene, cene+wene) );
  
  //3. inverse
  CpxNumTnsBlkd Y;
  iC( fdct3d_inverse(m,n,p, nbscales,nbdstz_coarse, C,W,Y) );
  tm1 = time(NULL);  iC( PetscPrintf(MPI_COMM_WORLD, "INVERSE %e\n", difftime(tm1,tm0)) );  tm0 = tm1;
  double yene = Y.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "Y energy %e\n", yene) );
  
  //4. compute difference
  for(int i=0; i<e; i++)	 for(int j=0; j<f; j++)		for(int k=0; k<g; k++) {
	 if(X.owners()(i,j,k)==mpirank) {
		iA(Y.owners()(i,j,k)==mpirank); 		//CHECK THEY HAVE THE SAME DISTRIBUTION
		CpxNumTns& Xblk = X.block(i,j,k);		CpxNumTns& Yblk = Y.block(i,j,k);
		for(int ioff=0; ioff<b; ioff++)		  for(int joff=0; joff<b; joff++)			 for(int koff=0; koff<b; koff++) {
		  Yblk(ioff,joff,koff) -= Xblk(ioff,joff,koff);
		}
	 }
  }
  double eene = Y.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "E energy %e\n", eene) );
  /*
  if(mpisize==1) {
	 for(int s=0; s<nbscales-1; s++) {
		int nw = C.owners()[s].size();
		for(int w=0; w<nw/2; w++) {
		  CpxNumTns& A = C.block(s,w);		  CpxNumTns& B = C.block(s,w+nw/2);
		  double maxerr = 0;
		  for(int i=0; i<A.m(); i++)			 for(int j=0; j<A.n(); j++)				for(int k=0; k<A.p(); k++)
			 maxerr = max(maxerr, abs(A(i,j,k)-conj(B(i,j,k))));
		  cerr<<s<<" "<<w<<" "<<maxerr<<endl;
		}
	 }
  }
  */

  PetscFinalize();
  return 0;
}
