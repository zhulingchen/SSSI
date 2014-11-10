function [vel,x,z]=dipmodel(dx,xmax,zmax,vhigh,vlow,dip,znot)
% dipmodel: build a simple 2D model of a dipping reflector
% 
% [vel,x,z]=dipmodel(dx,xmax,zmax,vhigh,vlow,dip,dnot)
%
% dx ... grid interval (distance between grid points in x and z)
% xmax ... maximum x coordinate (minimum is zero)
%  *********** default 2500 **********
% zmax ... maximum z coordinate (minimum is zero)
%  *********** default 1000 ************
% vhigh ... velocity below the dipping reflector
%  *********** default 4000 ************
% vlow ... velocity above the dipping reflector
%  *********** default 2000 ************
% dip ... dip in degrees of the reflector, positive is down to the right
%  *********** default = 20 degrees ************
% znot ... depth to dipping reflector at xmax/2
%  *********** default = zmax/4 **************
%
% vel ... velocity model matrix
% x ... x coordinate vector for vel
% z ... z coordinate vector for vel
%
% NOTE: the simplest way to plot vel is: plotimage(vel-mean(vel(:)),z,x)
%

if(nargin<5)
    vlow=2000;
end
if(nargin<4)
    vhigh=4000;
end
if(nargin<3)
    zmax=1000;
end
if(nargin<2)
    xmax=2500;
end
if(nargin<6)
    dip=10;
end
if(nargin<7)
    znot=zmax/3;
end

x=0:dx:xmax; % x coordinate vector
z=0:dx:zmax; % z coordinate vector

%initialize velocity matrix as a constant matrix full of vlow
vel=vlow*ones(length(z),length(x));

% define the dipping horizon as a three point polygon
dx2=dx/2;
xnot=xmax/2;
x1=x(1)-dx2;x2=x(end)+dx2;
z1=znot+(x1-xnot)*tand(dip);
z2=znot+(x2-xnot)*tand(dip);

xpoly=[x1 xnot  x2 x2 x1];zpoly=[z1 znot z2 zmax+dx2 zmax+dx2];

% install the basement in the velocity matrix
vel=afd_vmodel(dx,vel,vhigh,xpoly,zpoly);

