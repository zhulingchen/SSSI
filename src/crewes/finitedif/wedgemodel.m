function [vel,x,z]=wedgemodel(dx,xmax,zmax,vhigh,vlow)
% WEDGEMODEL: build a simple 2D model of a high-velocity wedge over an anticline
% 
% [vel,x,z]=wedgemodel(dx,xmax,zmax,vhigh,vlow)
%
% dx ... grid interval (distance between grid points in x and z)
% xmax ... maximum x coordinate (minimum is zero)
%  *********** default 2500 **********
% zmax ... maximum z coordinate (minimum is zero)
%  *********** default 1000 ************
% vhigh ... velocity in the wedge and the anticline
%  *********** default 4000 ************
% vlow ... velocity below the wedge and above the anticline
%  *********** default 2000 ************
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

xpinch=xmax*3/5; % wedge pinchout coordinates
zwedge=zmax/2; % wedge maximum depth
x=0:dx:xmax; % x coordinate vector
z=0:dx:zmax; % z coordinate vector

%initialize velocity matrix as a constant matrix full of vlow
vel=vlow*ones(length(z),length(x));

% define the wedge as a three point polygon
dx2=dx/2;
xpoly=[-dx2 xpinch -dx2];zpoly=[-1 -1 zwedge];

% install the wedge in the velocity matrix
vel=afd_vmodel(dx,vel,vhigh,xpoly,zpoly);

% define an anticline beneath the wedge

x0=xpinch/2;z0=zwedge+100; % x and z of the crest of the anticline
a=.0005; % a parameter that determines the steepness of the flanks
za=a*(x-x0).^2+z0; % model the anticline as a parabola

% build a polygon that models the anticline
%ind=near(za,zmax+dx);
ind=find(za<=zmax-20*dx);
xpoly=[x(ind(1))-dx x(ind(2:end)) x(ind(1))-dx ];zpoly=[za(ind) za(ind(end))];

%install the anticline in the velocity model
vel=afd_vmodel(dx,vel,vhigh,xpoly,zpoly);

% bottom layer
z1=max(zpoly)-dx;
xpoly=[-dx xmax+dx xmax+dx -dx];zpoly=[z1 z1 zmax+dx zmax+dx];

vel=afd_vmodel(dx,vel,vhigh,xpoly,zpoly);
