function dipspect(a,t,dx,f,vo,c,z,thrline,thapline,thaline)
% DIPSPECT: make dip spectral analysis plot for v(z)
%
% dipspect(a,t,dx,f,vo,c,z,thrline,thapline,thaline)
% dipspect(a,t,dx,f,vo,c,z)
%
% DIPSPECT makes a dip spectral analysis plot for v=vo+cz
% If z is not provided, then it is computed with 
% z=linspace(0,zmax,nz). Here, zmax is computed as the maximum
% depth for the supplied record length and nz=500. (If z is
% provided as a single number, then it is taken to be nz).
%
% a ... aperture
% t ... maximum recording time
% dx ... spatial sampling interval
% f ... frequency of interest
% vo ... initial velocity
% c ... accelerator ( v(z) = vo +c*z )
% z ... vector of depths for which limits are computed
% thrline ... color to plot the record length limit with
%    ********* default 'r' *********
% thapline ... color to plot the aperture limit with
%    ********** default 'c' ********
% thaline .... color to plot with alias limit with
%    ********* default 'g' ********
%
%
% G.F. Margrave, CREWES Project, 1997
%
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

if(nargin<10)
	thaline='g';
end
if(nargin<9)
	thapline='c';
end
if(nargin<8)
	thrline='r';
end

if(nargin<7)
	z=500;
end

if(length(z)==1)
	[thetar,z]=threc(t,vo,c,z);
else
	thetar=threc(t,vo,c,z);
end
ir=find(imag(thetar)==0);

thap=thaper(a,vo,c,z);
inda=find(imag(thap)==0);
thal = thalias(dx,f,vo,c,z);
indal=find(imag(thal)==0);

figure;

plot(z(ir),thetar(ir),thrline,z,thap,thapline,z(indal),thal(indal),thaline);
legend(['Record length limit, T=' num2str(t) ' sec'],...
	['Aperture limit, A=' int2str(a) ' lu'],...
	['Spatial aliasing limit, dx=' int2str(dx) ' lu, f=' ...
	int2str(f) ' Hz']);
	
title([' Dipspect chart for vo= ' int2str(vo) ' c= ' num2str(c)])
grid

