/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#include "cpxcrvletocr.hpp"

CpxCrvletOcr::CpxCrvletOcr(const char* name)
{
  strcpy(_name, name);
}

CpxCrvletOcr::~CpxCrvletOcr()
{
  //remove
  vector< vector<int> >& c = _nxs;
  for(int si=0; si<c.size(); si++)
	 for(int wi=0; wi<c[si].size(); wi++) {
		char tmp[100];		sprintf(tmp, "%s_%d_%d.dat", _name, si, wi);
		remove(tmp);
	 }
}

int CpxCrvletOcr::setup(vector< vector<int> > nxs, vector< vector<int> > nys, vector< vector<int> > nzs, int maxnb)
{
  _nxs = nxs;  _nys = nys;  _nzs = nzs;
  
  _maxnb = maxnb;
  _count = 0;
  _clock = 0;
  
  vector< vector<int> >& c = nxs;
  _blocks.resize(c.size());  _szvec.resize(c.size());  _tmvec.resize(c.size());
  for(int s=0; s<c.size(); s++) {
	 _blocks[s].resize(c[s].size());  _szvec[s].resize(c[s].size());  _tmvec[s].resize(c[s].size());
	 for(int w=0; w<c[s].size(); w++) {
		_blocks[s][w].resize(0,0,0); //empty initially
		_szvec[s][w] = _nxs[s][w]*_nys[s][w]*_nzs[s][w];
		_tmvec[s][w] = -1;
	 }
  }
  
  return 0;
}

CpxNumTns& CpxCrvletOcr::block(int s, int w)
{
  if(_blocks[s][w].m()==0) { //not in
	 if(_count==_maxnb) { //no space
		//kick someone out
		int low = _clock;
		int si = -1;		int wi = -1;
		vector< vector<int> >& c = _szvec;
		for(int s=0; s<c.size(); s++)
		  for(int w=0; w<c[s].size(); w++) {
			 if(_tmvec[s][w]!=-1)
				if(_tmvec[s][w]<low) {
				  low = _tmvec[s][w];
				  si = s;				  wi = w;
				}
		  }
		char tmp[100];		sprintf(tmp, "%s_%d_%d.dat", _name, si, wi);
		ofstream fout(tmp);
		fout.write((char*)(_blocks[si][wi].data()), _szvec[si][wi] * sizeof(cpx));
		fout.close();
		
		_blocks[si][wi].resize(0,0,0);
		_tmvec[si][wi] = -1;
		_count--;
	 }
	 //load in
	 _blocks[s][w].resize(_nxs[s][w], _nys[s][w], _nzs[s][w]);
	 char tmp[100];	 sprintf(tmp, "%s_%d_%d.dat", _name, s, w);
	 ifstream fin(tmp);
	 if(fin.good())
		fin.read((char*)(_blocks[s][w].data()), _szvec[s][w] * sizeof(cpx));
	 _tmvec[s][w] = _clock;
	 _count ++;
	 _clock ++;
  } else { //in
	 if(_clock-_tmvec[s][w] > _maxnb/2) {
		_tmvec[s][w] = _clock;		_clock ++;
	 }
  }
  return _blocks[s][w];
}

