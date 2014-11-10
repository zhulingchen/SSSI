function [xc,yc,radius] = twoptcircle(xp,yp,n)
% [xc,yc,radius] = twoptcircle(xp,yp,n)
% [xc,yc,radius] = twoptcircle(xp,yp)
%
% TWOPTCIRCLE takes two points, and makes a circle with point 1 as the
% center, and point 2 defining the radius.
%
% xp	= x-coordinates for point1 and point2, i.e. x(1) and x(2)
% yp	= y-coordinates for point1 and point2, i.e. y(1) and y(2)
% n	= number of points on the circle
% **************************** Default = 30 *******************************
%
% xc	= x-coordinates for points on circle
% yc	= y-coordinates for points on circle
%
%  T. N. BISHOP,  OCTOBER 1993,  CPTC CANADA
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
if nargin<3
	n=30;
end
radius = sqrt( (xp(2)-xp(1)).^2 + (yp(2)-yp(1)).^2 );
for i = 1:n
  t(i) = (i-1)*2.*pi/n;
end
t = [t,t(1)];
  xc=xp(1) + radius*cos(t);
  yc=yp(1) + radius*sin(t);
