function dmin = polydist(x,y,datum,icol)
% dmin = polydist(x,y,datum,icol)
% dmin = polydist(x,y,datum)
%
% POLYDIST is the minimum distance from a polynomial to a datum line
%
% x,y 	= coordinate vectors of polygon 
% datum = 4-vector, [x1 y1 x2 y2] definining 2 pts on datum
% icol 	= color for plotting, 'y','m','c','r','g', or 'b'
% ====================== Default = 'b' =======================
%
% T.N.Bishop  Oct.93
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
% get slope for datum
m = zeros(size(y));  %slope of datum
m = (datum(4)-datum(2))/(datum(3)-datum(1));
if(datum(3) == datum(1)) 
  m = 1/eps;
end
b = datum(2) - m*datum(1);
%find dist from each point in polygon to datum line
msq = m.^2;
denom = 1+msq;
xd=zeros(size(y));   %intercept of datum for each perpendicular
yd=zeros(size(y));   %y coord of intercept of datum 
d = zeros(size(y));  %dist.along perpendicular
%     perpendicular is perpendicular to the datum and it also
%     intersects the polygon point (x(i), y(i)).
xd =( x+m.*(y-b) )./denom;
yd =( x*m+msq*y+b )./denom;
d = sqrt((x-xd).^2 + (y-yd).^2);
dmin = min(d);
imin = find(d == dmin);
% plot segment (if icol is called by the function)
if(nargin > 3) 
  hold on
  xtmp=[x(imin) xd(imin)];
  ytmp=[y(imin) yd(imin)];
  plot(xtmp,ytmp,icol)
end
