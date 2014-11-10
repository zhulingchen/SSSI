function [x0,t0,dtdx,tray,r]=normray(xn,zn,dip,params,clr)
% NORMRAY: trace a normal ray to the surface
%
% [x0,t0,dtdx,tray,r]=normray(xn,zn,dip,params)
%
% 
% (xn,zn) ... coordinates of the normal incidence point
% dip ... reflector dip in degrees
% params(1) ... figure number to plot in (time)
%       *********** default = don't plot ********
% params(2) ... figure number to draw ray in depth
%       *********** default = don't plot ********
% params(3) ... time step
%       *********** default =.004 **********
% params(4) ... 1 for one-way times, 2 for two-way times
%       ************ default = 2 ***********
% params(5) ... maximum time allowed for raytracing
%       ************ determined automatically from velocity model ********
% clr ... color to draw with
%       ************ default = 'r' *********
%
% x0 ... lateral coordinate at which the ray arrives at z=0
% t0 ... time of the ray (to surface)
% dtdx ... time-dip (horizontal slowness) of the ray at surface (includes
%		the factor of 2 if two-way times are used)
% tray ... vector of times for the ray, at increments of params(3)
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

global RTV2 RTDLNVDX RTDLNVDZ RTDG RTVRMS RTVAVE RTT RTZ

if(nargin<5) clr='r'; end
if(nargin<4) params=nan*[1:5]; 
elseif (length(params)~=5)
	error('params must be specified with 5 entries')
end

if(isnan(params(1))) params(1) = 0; end
if(isnan(params(2))) params(2) = 0; end
if(isnan(params(3))) params(3) = .004; end
if(isnan(params(4))) params(4) = 2; end
if(isnan(params(5))) 
	%guess time for raytracing
    if(zn>max(RTZ))
        zn=max(RTZ);
    end
	t0=2*interp1(RTZ,RTT,zn)/(cos(pi*dip/180)+.01);
	params(5)=1.5*t0;
end


tfigno=params(1); zfigno=params(2); dt=params(3); timefact=params(4);

tmax=params(5); %maximum time to raytrace

%define ray
[nx,nz]=size(RTV2);
x=(0:nx-1)*RTDG;
z=(0:nz-1)*RTDG;
ix=near(x,xn);
iz=near(z,zn);
if(iz<1||iz>size(RTV2,1)||ix<1||ix>size(RTV2,2))
    x0=nan;
    t0=nan;
    dtdx=nan;
    tray=nan;
    r=nan;
    return;
end
vn=sqrt(RTV2(iz,ix));
p=sin(pi*dip/180)/vn;
q=-sqrt(1/vn^2 - p^2);
r0=[xn,zn,p,q];

	[tray,r]=shootraytosurf(dt,tmax,r0);
	zray=r(:,2);
	if(zray(end)>0)
		disp('ray reached max time before reaching surface')
		t0=tmax;
	else
		testz=0;
		if(zray(end)~=0)
			disp('interp needed');
		end
	end
 
tray=timefact*tray;
t0=max(tray);
x0=r(length(zray),1);
dtdx=timefact*r(length(zray),3);


%plot if requested
if(zfigno)
	figure(zfigno);
	hray=line(r(:,1),r(:,2),ones(size(r(:,1))),'color',clr);
	%compute a perpendicular
	n=length(zray);
	if(r(2,1)~=r(1,1))
		slope=(r(2,2)-r(1,2))/(r(2,1)-r(1,1));
	else
		slope=inf;
	end
	if(abs(slope)>100*eps)
		slopeperp=-1/slope;
	else
		slopeperp=inf;
	end
	x1=r(1,1)-RTDG;x2=r(1,1)+RTDG;
	z1=r(1,2)-RTDG*slopeperp;z2=r(1,2)+RTDG*slopeperp;
	hdip=line([x1 x2],[z1 z2],[1 1],'color',clr,'linewidth',2);
end
%plot if requested
if(tfigno)
	figure(tfigno);
	delx=6*RTDG;
	x1=x0-delx;x2=x0+delx;
	t1=t0-delx*dtdx;t2=t0+delx*dtdx;
	hdip=line([x1 x2],[t1 t2],[1 1],'color',clr,'linewidth',3);
end
	
		
