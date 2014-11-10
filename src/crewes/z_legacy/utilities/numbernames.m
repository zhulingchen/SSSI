function namesnum=numbernames(names)
% namesnum=numbernames(names)
%
% Function takes a matrix of names (string matrix) and searchs for
% duplicate names. Any duplicate names are suffixed with numbers to
% make them unique. The first occurance of any name is left without
% any suffix while the second has 02 and on up to 99. The returned 
% name matrix will have two extra columns which will contain the 
% numbers. (If no duplicate names are found, then the columns are
% stiil there with blanks.)
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
[n,c]=size(names);
namesnum=[names setstr(32*ones(n,2))];
for k=1:n-1
	testname=names(k,:);
	
	if(isempty(deblank(namesnum(k,c+1:c+2))))%see if its already flagged
	
		ind=strmatfind(names(k+1:n,:),testname);
	
		if(~isempty(ind))
			for kk=1:length(ind)
				snum=int2str(kk+1);
				l=length(snum);
				namesnum(k+ind(kk),c+3-l:c+2)=snum;
			end
		end
	end
end
