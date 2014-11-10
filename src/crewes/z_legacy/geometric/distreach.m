function [s,indicies]=distreach(x,y,xnot,ynot,r)
% [s,indicies]=distreach(x,y,xnot,ynot,r)
%
% This function searches the point set represented by x,y to find those points
% which lie within a radius r of the single point xnot,ynot. The distances of
% the found points are returned in s and their indicies in indicies. Thus, the
% number of found points is length(s) and they are identified as:
% x(indicies),y(indicies) is the set of points whose distances from xnot,ynot
% are given by s <= r.
% See also distcum distinc and distpoint
%
% x = vector containing the x coordinates of the data points
%
% y = vector containing the y coordinates of the data points
%
% xnot = scalar x coordinate of the single point
%
% ynot = scalar y coordinate of the single point
%
% r= radius of the search region (also called the reach)
%
% s = distance vector containing the distance between each point in
%     x and y and the point xnot,ynot
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
% do first sort
xmin = xnot-r;xmax=xnot+r;
ymin= ynot-r; ymax= ynot+r;
ind = find((x<=xmax)&(x>=xmin));
ind = find((y(ind)<=ymax)&(y(ind)>=ymin));
s=distpoint(x(ind),y(ind),xnot,ynot);
indicies = find( s<= r);
s=s(indicies);
indicies=ind(indicies); 
 
