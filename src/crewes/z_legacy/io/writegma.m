function flag=writegma(filename,logname,z,samps,units,kb,topnames,...
		ztops)
% flag=writegma(filename,logname,z,samps,units,kb,topnames,...
%				ztops)
%
% Write a log to disk in GMA format
% filename ... string containing the file name
% logname ... string containing the logname. If longer than 16
%	characters, it will be truncated.
% z ... vector of depths for the log
% samps ... vector of log samples. (z and samps must be the same
%		length)
% units ... Units flag. Use 'M' for metric and 'F' for imperial
% kb ... kelly bushing elevation
%	******* default = 0.0 ******
% topnames ... string matrix of top names. One name per row. 
%	Maximum name length is 30 (16 in LOGM)
%  ******** default is no tops *******
% ztops ... vector of formation top depths. Must be one per
%	topname.
%  ******** default is no tops *******
%
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.
 
% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by 
% its author (identified above) and the CREWES Project.  The CREWES 
% project may be contacted via email at:  crewesinfo@crewes.org
% 
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) Use of this SOFTWARE by any for-profit commercial organization is
%    expressly forbidden unless said organization is a CREWES Project
%    Sponsor.
%
% 2) A CREWES Project sponsor may use this SOFTWARE under the terms of the 
%    CREWES Project Sponsorship agreement.
%
% 3) A student or employee of a non-profit educational institution may 
%    use this SOFTWARE subject to the following terms and conditions:
%    - this SOFTWARE is for teaching or research purposes only.
%    - this SOFTWARE may be distributed to other students or researchers 
%      provided that these license terms are included.
%    - reselling the SOFTWARE, or including it or any portion of it, in any
%      software that will be resold is expressly forbidden.
%    - transfering the SOFTWARE in any form to a commercial firm or any 
%      other for-profit organization is expressly forbidden.
%
% END TERMS OF USE LICENSE
if(nargin<6)
	kb=0.0;
end
if(nargin<7)
	topnames=[];
	ztops=[];
end
%open the file
fid= fopen(filename,'wt');
if(fid==-1)
	flag=fid;
	error('Unable to open output file')
end
rec=blanks(80);
if(length(logname)>16) logname=logname(1:16); end
rec(1:length(logname))=logname;
zbeg=z(1);zend=z(length(z));dz=z(2)-z(1);
s=sprintf('%7.1f',zbeg);
if(length(s)>7) s=s(1:7); end
rec(29:28+length(s))=s;
s=sprintf('%7.1f',zend);
if(length(s)>7) s=s(1:7); end
rec(36:35+length(s))=s;
s=sprintf('%6.4f',dz);
if(length(s)>6) s=s(1:6); end
rec(47:46+length(s))=s;
s=sprintf('%8.0f',length(z));
rec(54:53+length(s))=s;
rec(70)=units;
%write first record
fprintf(fid,'%s\n',rec);
%loop to write out samples
rec=blanks(80);
inum=1;
nz=length(z);
while(inum < nz)
	rec(1:7)=sprintf('%7.1f',z(inum));
	nnums=min([10 nz-inum+1]);
	rec(8:7+nnums*7)=sprintf('%7.1f',samps(inum:inum+nnums-1));
	rec(78:80)=sprintf('%3d',nnums);
	fprintf(fid,'%s\n',rec);
	rec=blanks(80);
	inum=inum+nnums;
end
	
%write out tops
ntops=size(topnames,1);
if(ntops)
	rec=sprintf('#TOPS#%5d',ntops);
	fprintf(fid,'%s\n',rec);
	for k=1:ntops
		name=strunpad(topnames(k,:));
		if(length(name)>30) name=name(1:30); end
		rec=blanks(40);
		rec(1:length(name))=sprintf('%s',name);
		rec(33:39)=sprintf('%7.1f',ztops(k));
		fprintf(fid,'%s\n',rec);
	end
end
	%kb
	rec=sprintf('KB%7.1f',kb);
	fprintf(fid,'%s\n',rec);
