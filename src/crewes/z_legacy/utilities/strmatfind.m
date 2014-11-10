function ind=strmatfind(smat,s)
% ind=strmatfind(smat,s)
%
% STRMATFIND searches a string matrix, smat, for a row whos contents
% contain the same string as s. ind is returned as a vector of row
% numbers indicating which strings match.
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
ind=[];
% unpad s
s=strunpad(s);
s=deblank(s);
ii=find(s~=32);
if(isempty(ii))
	return;
end
s=s(ii(1):length(s));
[nrows,ncols]=size(smat);
%form smat into a big row vector with string boundary markers
smat=[smat '|'*ones(nrows,1)];
smat=smat';
smat=['|' smat(:)' '|'];
ii=find(smat==1);
if(~isempty(ii))
	smat(ii)=[];
end
ii=find(smat==0);
if(~isempty(ii))
	smat(ii)=[];
end
ibnd=find(smat=='|');
ii=findstr(smat,s);
%check each ii for full match
for k=1:length(ii)
	i=ii(k);
	ib=surround(ibnd,i);
	s1=smat(ibnd(ib)+1:ibnd(ib+1)-1);
	s1=deblank(s1);
	ik=find(s1~=32);
	s1=s1(ik(1):length(s1));
	%test
	if(strcmp(s1,s))
		ind=[ind ib];
	end
end
