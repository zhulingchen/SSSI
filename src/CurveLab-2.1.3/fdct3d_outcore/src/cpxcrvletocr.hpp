/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/


#ifndef _CPXCRVLETOCR_HPP_
#define _CPXCRVLETOCR_HPP_

#include "numtns.hpp"

class CpxCrvletOcr
{
protected:
  vector< vector<int> > _nxs, _nys, _nzs;
  char _name[100];
  int _maxnb;
  int _count;
  
  int _clock;
  vector< vector<CpxNumTns> > _blocks;
  vector< vector<int> > _szvec;
  vector< vector<int> > _tmvec;
  
public:
  CpxCrvletOcr(const char* name);  //  CpxCrvletOcr(const CpxCrvletOcr& D);
  ~CpxCrvletOcr();  //  CpxCrvletOcr& operator=(const CpxCrvletOcr& D);
  int setup(vector< vector<int> > nxs, vector< vector<int> > nys, vector< vector<int> > nzs, int ma);
  CpxNumTns& block(int s, int w);
  //access
  vector< vector<int> >& nxs() { return _nxs; }
  vector< vector<int> >& nys() { return _nys; }
  vector< vector<int> >& nzs() { return _nzs; }
};

#endif
