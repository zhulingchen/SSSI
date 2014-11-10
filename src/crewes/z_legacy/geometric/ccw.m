function flag = ccw(pt1,pt2,pt3)
% flag = ccw(pt1,pt2,pt3)
%
% CCW is an implementation of Sedgewick's algorithm (Algorithm's in C++,
% Robert Sedgewick, 1992, Addison-Wesley, p350) for determining the sense of
% rotation when traveling from pt1 to pt2 to pt3 (these are points in the
% x,y plane). If this travel results in rotation through a counter clockwise
% angle then +1 is returned which -1 is returned for a clockwise angle. A
% special case is when the three points are colinear. In this case, if pt1
% is between pts 2&3 then -1 is returned, if pt2 is between pts 1&3 then +1
% is returned, and if pt3 is between pts 1&2 then 0 is returned. 
% 
% Note that each point is a 2 element vector giving x first then y.
%
% G.F. Margrave January 1994
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
dx1=pt2(1)-pt1(1); dy1=pt2(2)-pt1(2);
dx2=pt3(1)-pt1(1); dy2=pt3(2)-pt1(2);
if( dx1*dy2 > dy1*dx2 ) flag=1; return; end
if( dx1*dy2 < dy1*dx2 ) flag=-1; return; end
if( (dx1*dx2<0) | (dy1*dy2<0) ) flag=-1; return; end
if( (dx1*dx1+dy1*dy1) < (dx2*dx2+dy2*dy2) ) flag=1; return; end
flag=0;
