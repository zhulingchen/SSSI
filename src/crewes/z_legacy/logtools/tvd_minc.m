function [x,y,z]=tvd_minc(zmd,azmi,incl)
% [x,y,z]=tvd_minc(zmd,azmi,incl)
%
% TVD_MINC uses the minimum curvature method to compute true 3-D coordinates
% of a wellbore path given information from a deviation survey. The method is
% documented in 'Directional Survey Calculation', by J.T. Craig and B.V. Randal
% found in 'Petroleum Engineer',March, 1976.
%
% zmd ... vector of measured depths
% azmi ... vector of azimuth angles in degrees
% incl ... vector of inclination angles in degrees
%
% G.F. Margrave
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
zmd=zmd(:);
azmi=azmi(:);
incl=incl(:);
%make sure we have a 0 zmd
pad=0;
if(zmd(1)~=0)
	zmd=[0;zmd];
	azmi=[0;azmi];
	incl=[0;incl];
	pad=1;
end
torad= pi/180.;
phi=torad*incl;
theta=torad*azmi;
npts=length(zmd);
k=2:npts;
kminus=1:npts-1;
cosd= cos(phi(k)-phi(kminus)) - sin(phi(kminus)).*sin(phi(k)).*...
		(1. - cos(theta(k)-theta(k-1)));
tand= sqrt( cosd.^(-2) -1.);
dl = atan(tand);
ind= abs(dl) > .00001;
fc=zeros(size(dl));
fc(~ind)=ones(sum(~ind),1);
fc(ind)= 2*tan(dl(ind)/2)./dl(ind);
sinphi= fc.*(sin(phi(k)).*sin(theta(k))+sin(phi(kminus)).*sin(theta(kminus)))/2;
cosphi= fc.*(sin(phi(k)).*cos(theta(k))+sin(phi(kminus)).*cos(theta(kminus)))/2;
x=zeros(size(zmd));
y=zeros(size(zmd));
z=zeros(size(zmd));
x(k)= cumsum( (zmd(k)-zmd(kminus)).*cosphi );
y(k)= cumsum( (zmd(k)-zmd(kminus)).*sinphi );
z(k)= cumsum( fc.*(zmd(k)-zmd(kminus)).*(cos(phi(k))+cos(phi(kminus)))/2.);
if(pad)
	x=x(k);
	y=y(k);
	z=z(k);
end
%flip x and y to conform with MINCOM
tmp=x;
x=y;
y=tmp;
