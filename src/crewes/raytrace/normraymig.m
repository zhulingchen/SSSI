function [tray,r]=normraymig(x0,t0,dtdx,params,clr)
% NORMRAYMIG: migrate a normal incidence ray
%
% [tr,r]=normraymig(x0,t0,dtdx,params,clr)
%
% NORMRAYMIG performs normal incidence raytrace migration for a velocity
% model that can vary with x and z. The model is built by defining a 
% matrix of velocity values and calling RAYVELMOD. SHOOTRAYVXY is called
% to trace the ray. If desired, then the ray may be plotted in any
% desired figure window.
%
% x0 ... surface x location of the ray measurement
% t0 ... the measured traveltime for the ray
% dtdx ... the measured time dip for the ray
% params(1) ... figure number for the ray to be plotted in
%    ********* default is no plotting *********
% params(2) ... time step size for raytracing (seconds)
%    ********* default = .004 seconds *********
% params(3) ... integer giving the number of grid increments
%		defining the width of the final dip bar on the plot.
%    ********* default = 4 ********
% PARAMS can be defaulted by providing nan for any entry.
% clr ... color to plot the rays with
%    ********* default = 'r' *************
% 
% tray ... vector of times from 0 to t0 sampled at params(2)
% r ... N by 4 matrix where N=length(tray). Each row is one time step. The
%	columns contain: r(:,1) x coordinate of the ray, r(:,2) z coordinate, 
%	r(:,3) horizontal slowness, r(:,4) vertical slowness.
%
% G.F. Margrave, CREWES, June 2000
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

global RTV2 RTDLNVDX RTDLNVDZ RTDG BILINEAR RTX

if(nargin<5) clr='r'; end
if(nargin<4) params=nan*[1:3]; end
if(isnan(params(1))) params(1)=0; end
if(isnan(params(2))) params(2)=.004; end
if(isnan(params(3))) params(3)=4; end
figno=params(1);
dt=params(2);
ndel=params(3);

x0=mean(x0);
t0=mean(t0);

%build x
%[nz,nx]=size(RTV2);
%x=RTDG*[0:nx-1];
x=RTX;

%setup
ix=near(x,x0);
v02=RTV2(1,ix);

% adjust step size to end exactly on tnot/2
tstep=0:dt:t0/2;
nsteps=length(tstep);
tstep=linspace(0,t0/2,nsteps);

p0=-.5*dtdx;q0=sqrt(1/v02-p0^2);
if(iscomplex(q0))
	disp('dtdx too steep, ray not physical');
	return;
end
r0=[x0,0,p0,q0];

[tray,r]=shootrayvxz(tstep,r0);

if(max(tray)<t0/2-1000*eps)
	disp('ray failed. Velocity model may be too small')
end

%plot if requested
if(figno)
	figure(figno);
	hray=line(r(:,1),r(:,2),ones(size(r(:,1))),'color',clr);
	%compute a perpendicular
	n=length(tray);
	if(r(n,1)~=r(n-1,1))
		slope=(r(n,2)-r(n-1,2))/(r(n,1)-r(n-1,1));
	else
		slope=inf;
	end
	if(abs(slope)>100*eps)
		slopeperp=-1/slope;
	else
		slopeperp=inf;
	end
	del=ndel*RTDG;
	x1=r(n,1)-del;x2=r(n,1)+del;
	z1=r(n,2)-del*slopeperp;z2=r(n,2)+del*slopeperp;
	hdip=line([x1 x2],[z1 z2],[1 1],'color',clr,'linewidth',2);
end
