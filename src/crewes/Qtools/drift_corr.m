function [sd,spd,tdr,tdp]=drift_corr(s,t,sp,z,Q,fs,f0)
% DRIFT_CORR ... apply drift corrections to a synthetic seismogram
%
% [sd,spd,tdr]=drift_corr(s,t,sp,z,Q,fs,f0) Use this if drift times are not known
% sd=drift_corr(s,t,tdr) use this mode if drift times are already known
%
% Sonic logging is done at roughly 12500Hz while seismic data typically has
% a dominant frequency below 50 Hz. Theory predicts that the velocities
% measured by the sonic tool will be systematically faster than those
% experienced by seismic waves. This frequency dependent velocity effect is
% the "dispersion" associated with Q attenuation. The implication is that
% the events on real data will be systematically delayed relative to those
% seen on a synthetic seismogram. This delay is called the drift. This
% effect is often addressed in practice by using a "check shot" survey and
% then adjusting the sonic log values. This function provides an
% alternative. If a synthetic has been created using the measured sonic
% values, then this function estimates and applies the drift delay to the
% synthetic.
%
% First mode
% s ... synthetic seismogram computed with measured sonic values
% t ... time coordinate for s (t(1) should be zero and corresponds to top
%       of log)
% sp ... P-wave sonic log used to compute s (not drift corrected)
%       Be sure this is a sonic and not the velocities from the sonic
% z  ... depth coordinate for sp
% Q  ... Q value for drift calculation. Can be either a scalar (single
%       value) or a vector the same size as sp.
% fs ... dominant frequency of seismogram (used for drift correction)
% f0 ... dominant frequency of sonic logging tool
%        ************ default f0 = 12500 Hz (12.5 kHz) ***********
%
% sd ... drift corrected synthetic
% spd ... drift corrected sonic log
% tdr ... drift time (same size as sp)
%
% Second mode
% s ... synthetic seismogram computed with measured sonic values
% t ... time coordinate for s (t(1) should be zero and corresponds to top
%       of log)
% tdr ... vector of drift times the same size as t. These are defined as
%       the traveltime at seismic frequency minus the traveltime at well
%       logging frequencies. They are typically positive numbers
%
% sd ... drift corrected synthetic
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

if( nargin>3)
    %check sonic
    smax=max(sp);
    if(smax>1000)
        error('It appears that your ''sonic'' is actually ''velocity''')
    end

    %check t(1)
    if(t(1)~=0)
        error('first time should be zero')
    end

    %check Q
    if(length(Q)>1)
        if(length(Q)~=length(z))
            error('Q must either be a scalar or a vector the same length as the sonic')
        end
    else
        Q=Q*ones(size(z));
    end

    if(nargin<7)
        f0=12500;
    end

    %adjust log so that first depth is 0. This ensures tdr is zero at start
    %of log
    z=z-z(1);

    %
    vwell=10.^6./sp;

    tdr=2*tdrift(Q,z,fs,vwell,f0);
    %fit with low-order poly with no constant term (this ensures a drift time
    %of zero at z=0
    % z0=z/mean(z);
    % M=zeros(length(z0),m);
    % for k=1:m
    %   M(:,k)=z0.^k;
    % end
    % p=M\tdr;
    % tdp=M*p;

    %compute a time-depth curve from the sonic
    [tz,zt]=sonic2tz(sp,z,-100,0);

    %convert the drift curve from depth to time
    t_td=interp1(zt,tz,z);%times for tdr

    %interpolate values from tdr at sonic times
    td_s=interp1(t_td,tdr,t);

    %now, drift correct the sonic.
    vs=(2*z.*vwell)./(2*z+tdr.*vwell);
    %check for nan in first place
    if(isnan(vs(1)))
        vs(1)=vwell(1);
    end
    if(sum(isnan(vs)))
        error('Drift correction logic failure')
    end
    spd= 10^6 ./vs;%drift corrected sonic values
else
    td_s=sp;
end

%ok, so we want a new synthetic where the sample at time t comes from the
%original synthetic at the earlier time t-td_s
t_int=t-td_s;%interpolation sites
ind=find(t_int<0);
if(~isempty(ind))
    t_int(ind)=0;
end
%sinc function interpolation
sd=sinci(s,t,t_int);