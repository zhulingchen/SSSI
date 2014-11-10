function ind=surround(x,xtest)
% SURROUND: analyze how a vector surrounds some test points
%
% ind=surround(x,xtest)
%
% SURROUND returns a vector of indicies indicating how the vector
% x surrounds xtest which must be a scalar. If 
% isempty(ind) then xtest lies outside the range of x. Otherwise,
% ind will be the index of a point in x just greater (or less) than 
% xtest. Thus the following will be true:
%	x(ind) <= xtest < x(ind+1)
%		or
%	x(ind) >= xtest > x(ind+1)
%
% So, for if xtest is an interior point for the vector x,
% ind and ind+1 select those points in x which surround (or bracket)
% xtest. Note that x need not be monotonic. If the xtest is surrounded
% more than once by x, then ind will be a vector.
% example: x=1:10;
% >>surround(x,-1)
%    returns []
% >>surround(x,3)
%   returns 3
% now let x=[1:10 9:-1:1]
% >>surround(x,3)
%    returns 3 17
% >>surround(x,pi)
%   returns 3 16
%
% G.F. Margrave
% Jan 1994, revised Jan 95
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

	
	n=length(x);
	x1=x(1:n-1);
	x2=x(2:n);
	
	ind=find( (x1<=xtest & x2>xtest ) | ...
			(x1>=xtest & x2 < xtest ) );

