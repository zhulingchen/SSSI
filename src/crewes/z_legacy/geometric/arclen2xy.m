function [xnot,ynot]=arclen2xy(x,y,a)
% [xnot,ynot]=arclen2xy(x,y,a)
%
% Given a piecewise linear curve whose nodes are described by the vectors
% x & y, ARCLEN2XY computes the (x,y) coordinates of the points whose
% arclength (or inline distance) from (x(1),y(1)) is given by the vector a.
% If the arclength is too large (or negative) then NaN is returned.
%
% G.F. Margrave December 1993
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
anodes=xy2arclen(x,y);
anodes(length(anodes))=anodes(length(anodes))+5*eps;
n=length(a);
xnot=nan*ones(size(a));
ynot=nan*ones(size(a));
for k=1:n
	ind=surround(anodes,a(k));
	if(~isempty(ind))
		ind=ind(1);% make sure we take the first one
		if( ind>0 )
			factor=(a(k)-anodes(ind))/(anodes(ind+1)-anodes(ind));
			xnot(k)=x(ind)+factor*(x(ind+1)-x(ind));
			ynot(k)=y(ind)+factor*(y(ind+1)-y(ind));
		end
	end
	
end
