function rayvelmod(v,dg,z,x)
% RAYVELMOD: establish a velocity model for vxz raytracing
%
% function rayvelmod(v,dg,z,x)
%
% Establish velocity model for raytracing
% v ... 2D matrix of velocities
% dg ... grid size for v. (i.e. both dx and dz are dg)
% z ... depth coordinate vector for v
% ******* default dg*(0:size(v,1)-1)' ***********
% x ... lateral coordinate vector for v
% ******* default dg*(0:size(v,2)-1) ***********
%
% This function stores the velocity model as a global. If
% you wish to save space, you can delete your model after
% calling RAYVELMOD.
%
% Globals established
% RTV2 ... matrix of velocity squared
% RTDLNVDX ... matrix of d(ln(v))/dx
% RTDLNVDZ ... martix of d(ln(v))/dz
% RTDG ... grid spaceing (x&z) in physical units
% RTVAVE ... a single average velocity function (lateral mean)
% RTVRMS ... a single rms velocity function (lateral mean)
% RTT ... one-way time coordinate for RTVAVE and RTVRMS
% RTZ ... z coordinate for velocity model
% RTX ... x coordinate for velocity model
%
% global RTV2 RTDLNVDX RTDLNVDZ RTDG RTVRMS RTVAVE RTT RTZ
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

global RTV2 RTDLNVDX RTDLNVDZ RTDG RTVRMS RTVAVE RTT RTZ RTX

RTDG=dg;

if(nargin<3)
   z=((0:size(v,1)-1)')*dg;
end
if(nargin<4)
   x=(0:size(v,2)-1)*dg;
end

z=z(:);
x=x(:)';
   
%pad velocity model with four cells all around
npad=4;
nz=size(v,1);nx=size(v,2);
vp=zeros(nz+2*npad,nx+2*npad);
vp(npad+1:nz+npad,npad+1:nx+npad)=v;
%top
vp(1:npad,:)=ones(npad,1)*...
		[v(1,1)*ones(1,npad) v(1,1:nx) v(1,nx)*ones(1,npad)];
%bottom
vp(nz+npad+1:nz+2*npad,:)=ones(npad,1)*...
		[v(nz,1)*ones(1,npad) v(nz,1:nx) v(nz,nx)*ones(1,npad)];
%left
vp(npad+1:nz+npad,1:npad)=v(:,1)*ones(1,npad);
%right
vp(npad+1:nz+npad,nx+npad+1:nx+2*npad)=v(:,nx)*ones(1,npad);

RTV2=vp.^2;

[RTDLNVDX,RTDLNVDZ]=gradient(log(vp),dg,dg);

vm=mean(vp')';
%nz=length(vm);
RTZ=[z(1)-dg*((npad:-1:1)');z;z(nz)+dg*((1:npad)')];
RTX=[x(1)-dg*(npad:-1:1) x x(nx)+dg*(1:npad)];

RTT=vint2t(vm,RTZ);
RTVAVE=vint2vave(vm,RTT);
RTVRMS=vint2vrms(vm,RTT);

disp('velocity model successfully converted to global variables');
