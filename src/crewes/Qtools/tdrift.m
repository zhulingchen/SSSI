function tdr=tdrift(Q,z,f1,v0,f0)
% TDRIFT: calculate the drift time
%
% tdr=tdrift(Q,z,f1,v0,f0)
%
% Sonic logging is done at roughly f0=12500Hz while seismic data typically has
% a dominant frequency (called f1) below 50 Hz. Theory predicts that the velocities
% measured by the sonic tool will be systematically faster than those
% experienced by seismic waves. This frequency dependent velocity effect is
% the "dispersion" associated with Q attenuation. The drift time is the
% traveltime difference using the velocities at f0 compared to f1. The
% formula used here is based on the Aki-Richards equation (5.81)
% v(f1)=v(f0)*(1+(1/(pi*Q))*log(f0/f1);
% which prescibes the variation of velocity with frequency for a constant Q
% medium. If t0 are the traveltimes at f0 and t1 are the traveltimes at f1,
% then the drift time is tdr=t0-t1. This will be positive if f0>f1.
%
% Q ... Q value (may be a vector)
% z ... distance traveled (may be a vector)
% f1 ... frequency of interest
% v0 ... velocity measured by logging tool at frequency f0
% f0 ... logging frequency (typically f1<<f0)
%   ********* default f0=12500 Hz ***********
% 
% tdr ... vector of drift times the same size as x
% Note: tdr is a one-way drift time. For application to surface reflction
% data this should be doubled.
%
%
% by G.F. Margrave, 2013
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

if(nargin<5)
    f0=12500;
end
if(length(Q)>1)
    if(size(Q)~=size(z));
        error('Q and Z must be the same size')
    end
end
%correct the velocities
v1=v0.*(1+(1./(pi*Q))*log(f1/f0));
%compute times at f0
t0=vint2t(v0,z);
%compute times at f1
t1=vint2t(v1,z);
tdr=t1-t0;%will be a positive quantity if f0>f1
