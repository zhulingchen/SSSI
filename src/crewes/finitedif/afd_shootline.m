function [shots,t,xshots,xrecs,shotnames]=afd_shootline(dx,vel,dt,dtstep,w,tw,tmax,nshots)
% AFD_SHOOTLINE ... shoot a 2D seismic line over a numerical model
%
% [shots,t,xshots,xrec]=afd_shootline(dx,vel,dt,dtstep,w,tw,tmax,nshots)
%
% dx ... spatial grid size of the velocity model (grid cell must be square)
% vel ... velocity matrix
% dt ... desire output time sample rate. This should usually be .004. It is
%        usually a waste of time to try for more out of a finite difference
%        simulation.
% dtstep ... time step size for the finite difference simulation. For
%        stability, it must be less than sqrt(3/8)*dx/max(vel(:)).
% w ... desired wavelet. Can be causal or non causal. Create this with
%       something like, wavemin, ricker, ormsby, wavez, or similar. It must be
%       sampled at the rate dt (not dtstep).
% tw ... time coordinate vector for w
% tmax ... maximum desired record time. If in doubt, run with a single shot
%       record first and check the results.
% nshots ... desired number of shot records. These will be placed at
%       regular intervals across the model. The actual locations are determined
%       by the funtion nominal_line_geom. Unless you want to wait a long
%       time, make this something like 20 or less.
%
% shots ... cell array of output shot gathers
% t ... time coordinate for the output traces
% xshots ... vector of shot locations
% xrecs ... cell array of receiver locations. Note that all shots will have
%           the same receivers positions. That is we simulate laying our a line of
%           receivers and shooting through it.
%
%
% G.F. Margrave, CREWES Project, August 2013
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

%we use 4th order laplacian. Test for stability
dt_test=sqrt(3/8)*dx/max(vel(:));
if(dt_test<dtstep)
    error(['dtstep too large, for stability it must be less than ' num2str(dt_test)]);
end
if((tw(2)-tw(2)-dt)>1000*eps)
    error('Wavelet must be sampled at the rate dt')
end
nx=size(vel,2);
x=(0:nx-1)*dx;
xmax=max(x);

%determine geometry
offedge=.05*xmax;
spreadlength=x(end)-x(1);
[xshots,xrecs]=nominal_line_geom(x,nshots,spreadlength,offedge);

shots=cell(1,nshots);
shotnames=shots;

for k=1:nshots
    snap1=zeros(size(vel));
    snap2=snap1;
    snap2(1,round(xshots(k)/dx)+1)=1;%unit impulse at the shot position   
    [shots{k},shot,t]=afd_shotrec(dx,dtstep,dt,tmax,vel,snap1,snap2,xrecs{k},zeros(size(xrecs{k})),w,tw,2);
    shotnames{k}=['Shot ' int2str(k)];
    disp(['Completed shot ' int2str(k) ' of ' int2str(nshots)]);
end
