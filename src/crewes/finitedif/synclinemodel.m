function [vel,x,z]=synclinemodel(dx,xmax,zmax,vhigh,vlow,zsyncline,zfocal,radius)
% synclineMODEL : build a model representing a syncline in a stratigraphic sequence
%
% [vel,x,z]=synclinemodel(dx,xmax,zmax,vhigh,vlow,zsyncline,zfocal,radius)
%
% This function builds a velocity matrix representing a syncline beneath a
% homogeneous layer. The syncline is simply modelled as a semi-circle
% connected to a horizontal interface
%  
%                       x <== focal point of the syncline (zfocal)
%
% ------------*                 *--------------- <== Horizontal interface (zsyncline)
%              *               *<== Syncline
%                 *         *
%                    *****
% By choosing the focal point (center) of the circle below the surface
% (i.e. at positive depth), a burried focus with a re=verse-time branch
% will result.  Choosing the focus above the surface (negative depth) will
% show only a temporal syncline.
%
% dx ... grid interval (distance between grid points in x and z)
% xmax ... maximum x coordinate (minimum is zero)
%  *********** default 2500 **********
% zmax ... maximum z coordinate (minimum is zero)
%  *********** default 1000 ************
% vlow ... Velocity above the syncline
%  *********** default 2000 ************
% vhigh ... Velocity below the syncline
%  *********** default 4000 ************
% zsyncline ... depth to the horizontal interface surrounding the syncline
%  *********** default zmax/3 *********
% zfocal ... depth to the focal point of the syncline. This should not be
%       greater than zsyncline.
%  *********** default zmax/10 **********
% radius ... radius of the semi-circle defining the syncline
%  *********** default zmax/5 ***********
%
%NOTE: zfocal+radius must be greater than zsyncline or there will be no
%syncline.
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
    zsyncline=zmax/3;
end
if(nargin<7)
    radius=zmax/2;
end
if(nargin<8)
    zfocal=zmax/10;
end

% x and z coordinate vector
x=0:dx:xmax;
z=0:dx:zmax; 

%flood with vhigh
vel=ones(length(z),length(x))*vhigh;

%install the layer
x1=min(1)-dx;
x2=max(x)+dx;
z1=-dx;
z2=zsyncline;
xpoly=[x1 x2 x2 x1];
zpoly=[z1 z1 z2 z2];
vel=afd_vmodel(dx,vel,vlow,xpoly,zpoly);

%install syncline. Define a polygon roughly equivalent to the circle
npoints=500;
xf=xmax/2;%x coordinate of the focal point
xpoly=linspace(xf-radius,xf+radius,npoints);
zpoly=sqrt(radius^2-(xf-xpoly).^2)+zfocal;
vel=afd_vmodel(dx,vel,vlow,xpoly,zpoly);
