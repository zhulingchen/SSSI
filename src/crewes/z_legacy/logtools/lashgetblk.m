function [block,ibeg,iend]=lashgetblk(lash,bname)
% [block,ibeg,iend]=lashgetblk(lash,bname)
% 
% LASHGETBLK interrogates a LAS header to find a logical data
% 'block'. A block is defined as beginning with a line: 
% ~bname
% where bname is a string, and ending at the beginning of the next 
% block or the EOF. Note that the block is considered matched if
% the line in the LAS header matches ~bname up to the length of bname.
% Thus a line in an LAS header like 
% ~logedit parameters
% will match with bname of ~log, or ~logedit, etc... The returned block
% is the first matching one. The match is not case sensitive.
% lash ... LAS header in a string matrix
% bname ... string containing the block name
% block ... string matrix containing the block including the ~bname line.
% ibeg ... the line number in lash of the beginning of the block
% iend ... the line number in lash of the end of the block
%
% G.F. Margrave, Department of Geology and Geophysics,
%	University of Calgary, 1996
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
[nlines,m]=size(lash);
k=1;
ibeg=0;
iend=0;
target=['~' lower(bname)];
ln=length(target);
while k<=nlines
	% check for ~
	if(lash(k,1)=='~')
		str=lash(k,:);
		
		if(strcmp(lower(str(1:ln)),target) )
			ibeg=k;
		else
			if(ibeg)
				iend=k-1;
				break;
			end
		end
	end
	
	k=k+1;
end
if(iend & ibeg)
	block=lash(ibeg:iend,:);
else
	block=[];
end
	
