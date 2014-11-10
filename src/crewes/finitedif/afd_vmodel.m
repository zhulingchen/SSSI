function vmodout=afd_vmodel(dx,vmodin,vel,xpoly,zpoly);
% AFD_VMODEL ... makes simple polygonal velocity models
%
% vmodout=afd_vmodel(dx,vmodin,vel,xpoly,zpoly);
%
% AFD_VMODEL will superimpose a polygon with a different velocity onto
% the background velocity model.  The background model may be homogenous,
% layered, or as complicated as desired.  The program will return the
% velocity model with the polygon superimposed.
% HINT:  Plot your initial velocity model using imagesc with proper x and z 
% coordinates and then use ginput to pick the points of the polygon with your
% mouse and return the coordinates of the points.   
%
% dx = the bin spacing for both horizontal and vertical (in consistent units)
% vmodin = the 'background" velocity matrix in consistent units
%         the upper left corner is (0,0).
% vel = the velocity within the polygon in consistent units (scalar)
% xpoly = a vector of the x coordinates of the polygon
%       = can be entered in consistent units or bin numbers
% zpoly = a vector of the z coordinates of the polygon
%       = can be entered in consisent units or in bin numbers
%         NOTE:  the program will trace out the polygon
%         in the order of the coordinates entered - the 
%         order of your coordinates can effect the shape
%         of the polygon.
%
% vmodout = the velocity matrix with the polygon superimposed
%
% by Carrie Youzwishen, February 1999
% completely rewritten by G.F. Margrave June 2000
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

if(size(xpoly)~=size(zpoly))
	error('xpoly must be the same size as zpoly')
end
if(length(vel)~=1)
	error('vel must be a scalar')
end

vmodout=vmodin;
[nz,nx]=size(vmodin);
x=(0:nx-1)*dx;
z=((0:nz-1)*dx);
%cpts=ones(nz,1)*x+i*(z')*ones(1,nx);
%cpoly=xpoly+i*zpoly;

%in=insidepoly(cpts(:),cpoly.');
z=z';
xx=ones(size(z))*x;
zz=z*ones(size(x));
%in=inpolygon(xx,zz,xpoly,zpoly);
in=inpolygon(xx,zz,xpoly,zpoly);
%vmodout(on)=vel;
vmodout(in)=vel;