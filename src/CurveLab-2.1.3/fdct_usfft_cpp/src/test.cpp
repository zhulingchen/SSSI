/*
    Copyright (C) 2004 Caltech
   Written by Lexing Ying
*/

#include "fdct_usfft.hpp"
#include "fdct_usfft_inline.hpp"

using namespace std;
using namespace fdct_usfft_ns;

int optionsCreate(const char* optfile, map<string,string>& options)
{
  options.clear();
  ifstream fin(optfile); assert(fin.good());
  string name;  fin>>name;
  while(fin.good()) {
	 char cont[100];	 fin.getline(cont, 99);
	 options[name] = string(cont);
	 fin>>name;
  }
  fin.close();
  return 0;
}

int main(int argc, char** argv)
{
  clock_t ck0, ck1;
  srand48( (long)time(NULL) );

  assert(argc==2);
  //get options
  map<string, string> opts;  optionsCreate(argv[1], opts);
  
  //get input data
  map<string,string>::iterator mi;
  
  int m;
  mi = opts.find("-m"); assert(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>m; }
  int n;
  mi = opts.find("-n"); assert(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>n; }
  int nbscales;
  mi = opts.find("-nbscales"); assert(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>nbscales; }
  int nbangles_coarse;
  mi = opts.find("-nbangles_coarse"); assert(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>nbangles_coarse; }
  int ac;
  mi = opts.find("-ac"); assert(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>ac; }
  
  //generate mat
  srand48( (long)time(NULL) );
  CpxNumMat x(m,n);
  for(int i=0; i<m; i++)
	 for(int j=0; j<n; j++)
		x(i,j) = cpx(drand48(), 0);
  
  ck0 = clock();
  
  //fdct_usfft_
  vector< vector<CpxNumMat> > c; //c = Ax
  fdct_usfft(m, n, nbscales, nbangles_coarse, ac, x, c);
  ck1 = clock();  cout<<"FDCT_USFFT_  takes "<<double(ck1-ck0)/CLOCKS_PER_SEC<<" seconds "<<endl;  ck0 = ck1;
  
  //afdct_usfft_
  CpxNumMat s(m,n); //s = A*Ax
  afdct_usfft(m, n, nbscales, nbangles_coarse, ac, c, s);
  ck1 = clock();  cout<<"AFDCT_USFFT_ takes "<<double(ck1-ck0)/CLOCKS_PER_SEC<<" seconds "<<endl;  ck0 = ck1;
  
  //ifdct_usfft_
  CpxNumMat y(m,n);
  ifdct_usfft(m, n, nbscales, nbangles_coarse, ac, c, y);
  ck1 = clock();  cout<<"IFDCT_USFFT_ takes "<<double(ck1-ck0)/CLOCKS_PER_SEC<<" seconds "<<endl;  ck0 = ck1;
  
  CpxNumMat e(m,n);
  for(int i=0; i<m; i++)
	 for(int j=0; j<n; j++)
		e(i,j) = x(i,j) - y(i,j);
  cerr<<"accuracy of inversion "<<sqrt(energy(e)/(m*n))<<endl;
  
  //sum0, (x, x)
  cpx a0 = energy(x);
  //sum1, (s, x)
  cpx a1(0,0);  //a1.re = 0;  a1.im = 0;
  for(int i=0; i<m; i++)
	 for(int j=0; j<n; j++) {
		cpx sv = s(i,j);		cpx xv = x(i,j);
		a1 += conj(sv)*xv;
	 }
  //sum2, (c, c)
  cpx a2(0,0);  //a2.re = 0;  a2.im = 0;
  for(int sc=0; sc<c.size(); sc++)
	 for(int w=0; w<c[sc].size(); w++)
		a2 += energy(c[sc][w]);
  //cerr<<"x,x "<<a0<<endl;  //cerr<<"s,x "<<a1<<endl;  //cerr<<"c,c "<<a2<<endl;
  
  return 0;
}

