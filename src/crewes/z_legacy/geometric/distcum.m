function s=distcum(x,y)
% s=distcum(x,y)
% s=distcum(x)
%
% If the first form is used, x and y represent the 2-D coordinates of a set of
% points. If the second form is used, x represents an incremental distance 
% vector computed from the point coordinates by distinc. The return value, s,
% is the cumulative distance for the set of points. That is, s(j) represents
% the distance measured along the point set from the jth point to the first
% point.
% See also distpoint and distreach
%
% x = vector containing the incremental distances as computed by distinc
%
% s = cumulative distance vector containing the distance between each point
%      and the first point. 
%     Note: length(s) == length(x). s(1) == 0, s(2) == ( distance from
%     point 2 to point 1), s(j) = ( distance from point j to point 1 )
%
% by G.F. Margrave, March 1993
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
% 
if nargin > 1
	d = distinc(x,y);
else
	d=x;
end
s=zeros(size(d));
for j=2:length(d)
	s(j) = s(j-1)+d(j);
end
 
 
