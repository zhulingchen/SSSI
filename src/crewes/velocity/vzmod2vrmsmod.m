function [vrms,t]=vzmod2vrmsmod(vel,z,dt,tmax,flag)
%  VZMOD2VRMSMOD: Compute Vrms(x,t) from Vint(x,z)
%
%  [vrms,t]=vzmod2vrmsmod(vel,z,dt,tmax)
%
%  VZMOD2VRMSMOD converts an interval velocity model in depth (such as the
%  models required for finite difference modelling) into an rms velocity
%  model in time (such as is required for Kirchhoff time migration)
%
%  vel..........is the input velocity model in depth. Each row is a
%               constant depth.
%  z........depth coordinate for vel, length(z) = size(vel,1) 
%  dt.... desired time sample rate
%  tmax ... maximum two-way time desired
%  flag ... 1 ... extend the rms velocity to tmax by constant extrapolation
%           2 ... extend the interval velocity model to tmax by constant
%           extrapolation
%  ********** default flag = 1 ************
% NOTE: flag only matters if the maximum 2-way traveltime in the model is
% less than tmax. In that case, flag=1 extends the final vrms to tmax while
% flag=2 extends the final interval velocity to tmax and then computes
% vrms. The second option may seem more physical but can result in a very
% fast bottom layer.
%
%  vrms....is the output RMS velocity matrix in time
%  t....... output two-way time coordinate
%
%
%  Zoron Rodriguez, November 2006
%  G.F. Margrave, March 2013
%
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.

% BEGIN TERMS OF USE LICENSE
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

if(nargin<4)
    error('not enough input arguments')
end
if(nargin<5)
    flag=1;
end
if(length(z)~=size(vel,1))
    error('vel and z sizes are not compatible')
end
nx=size(vel,2);
t=(0:dt:tmax);%two way time for output
vrms=zeros(length(t),nx);
for k=1:nx
   tv=2*vint2t(vel(:,k),z);%two way time at kth location
   if(tv(end)<tmax && flag==2)%extend the last interval velocity if needed
       tv(end)=tmax;
   end
   vrms(:,k)=vint2vrms(vel(:,k),tv,t);
   if(rem(k,100)==0)
       disp(['finished location ' int2str(k) ' of ' int2str(nx)])
   end
end




