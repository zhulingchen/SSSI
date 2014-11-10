function areaavg = polyarea(xpoly,ypoly,n,icolsym)
% areaavg = polyarea(xpoly,ypoly,n,icolsym)
% areaavg = polyarea(xpoly,ypoly,n)
%
% POLYAREA computes the area of the polygon defined by
% arrays xpoly and ypoly (the x and y coordinates)
% n is the number of grid cells from xmin to xmax
% or nc, as refered to below.  nr=ymax-ymin is set to nc
%
% icolsym = 'color-symbol' to plot fill with, 'b+','c*', etc.
% ====================== Default = 'b+' ====================
%
%  Tom Bishop, Oct.93
%
%     TIMING
%     note that time is proportional to nc*nc
%     also proportional to length(xpoly)
%     takes 22 sec.on Sparc10 for nc=20, length=80
%     takes 11 sec.on Sparc10 for nc=20, length=40
%     ACCURACY
%     for nc=10, std is about 10%
%     for nc=20, std is about 3%
%     for nc=40, std is about 1%
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
if(nargin < 4)
  icolsym='c+';
end
if(nargin < 3)
  nc=20;
else
  nc=n;
end
nr=nc;    %make vertical grid cells same as horizontal
dx=(max(xpoly)-min(xpoly))/(nc-1);
dy=(max(ypoly)-min(ypoly))/(nr-1);
xx(1)=min(xpoly)-dx*1.5;
yy(1)=min(ypoly)-dy*1.5;
xx(2)=max(xpoly)+dx*.5;
yy(2)=max(ypoly)+dy*.5;
ia = 0;
dx=0;
dy=0;
%    note that we compute the area three times, original
%    grid, then shift grid by +dx/2, then shift by +dy/2
for loop1=1:2
  xmin=xx(1)+(loop1-1)*dx/2;
  xmax=xx(2)+(loop1-1)*dx/2;
  for loop2=loop1:2
    ymin=yy(1)+(loop2-1)*dy/2;
    ymax=yy(2)+(loop2-1)*dy/2;
    x = linspace(xmin,xmax,nc);
    y = linspace(ymin,ymax,nr)';
    xg = ones(size(y))*x;
    yg = y*ones(size(x));
    ind = inside(xg,yg,xpoly,ypoly);
    ninside=sum(ind);
    dx=(xg(1,2)-xg(1,1));
    dy=(yg(2,1)-yg(1,1));
    ia=ia+1;
    sumarea(ia)=ninside*dx*dy;
  end
end
areaavg=mean(sumarea);
if(nargin > 3)
  hold on
  plot(xg(ind),yg(ind),icolsym);
end
