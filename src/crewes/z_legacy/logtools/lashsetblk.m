function lashout=lashsetblk(lash,bname,block)
% lashout=lashsetblk(lash,bname,block)
%
% LASHSETBLK is the opposite side of the coin to LASHGETBLK. It searches
% an LAS header to find a logical data block named bname and replaces
% it with the string matrix called block. (If isempty(bname) then a deletion
% will occur. If the block does not already exist in the header then it
% it is inserted just before the ~A block (which begins the log traces). The
% ~A line MUST be the last line in the header or the header is invalid.
% lash ... string matrix containing the input las header
% bname ... string containing the block name
% block ... string matrix containing the block to be set (inserted)
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
%test ~A
[nlines,m]=size(lash);
if( ~strcmp(lower(lash(nlines,1:2)),'~a') )
	error('Invalid LAS header. The last line must begin with ~A')
end
if(isempty(block))
	block='';
end
%search for the block
[oldblock,ibeg,iend]=lashgetblk(lash,bname);
if(isempty(oldblock))
	top=lash(1:nlines-1,:);
	bot=lash(nlines,:);
else
	top=lash(1:ibeg-1,:);
	bot=lash(iend+1:nlines,:);
end
lashout=str2mat(top,block,bot);
	
