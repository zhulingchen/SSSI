function [ztrc,z]=time2depth(trc,t,tzcurve,dz)
% TIME2DEPTH: Convert a trace from time to depth by a simple stretch
%
% [ztrc,z]=time2depth(trc,t,tzcurve,dz)
%
% TIME2DEPTH converts a single trace from time to depth. The conversion is
% specified by a time-depth curve that is stored in a two-column matrix.
% The first column is a list of times and the second is a list of 
% corresponding depths. The times can be either one-way or two-way and
% need only be consistent with the time-coordinate vector of the input.
% Between points in the time-depth curve, depths are interpolated linearly.
% Note that this function does not apply an antialias filter. It is
% up to the user to ensure that dz is sufficiently small to preclude
% aliasing.
% 
% trc ... the trace in time
% t ... time coordinate vector for trc
% tzcurve ... an n-by-2 matrix giving the time-depth curve.
%	n is the number of points on the curve and is arbitrary.
%	(more points is usually more accurate.) The first column is
%	time and the second column is the depths that correspond to
%	the times in the first column. The first time should be zero
%	and the last should be greater than or equal to max(t). This can be
%	created with sonic2tz.
% dz ... depth sample size. 
%
% NOTE: to avoid aliasing pick dz<vmin*dt/n where n is 1 for one-way time
% and 2 for 2way time and vmin is the slowest velocity in the model.
%
% G.F. Margrave, CREWES, Nov, 2000
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


%check tz curve limits
tz=tzcurve(:,1);zt=tzcurve(:,2);
if(min(tz)~=0)
	error('tz curve must start at zero time');
end
if(max(tz)<max(t))
	error('tz curve must extend to times greater than max(t)');
end

%make sure depths are monotonic
ztest=diff(zt);
ind=find(ztest<=0, 1);
if(~isempty(ind))
	error('depths on tzcurve must increase monotonically')
end

%determine depth range
z1=pwlint(tz,zt,t(1));
z2=pwlint(tz,zt,max(t));

z=((0:round((z2-z1)/dz))*dz)';

%ztrc=zeros(size(z));

%interpolation sites
tint=pwlint(zt,tz,z);

%sinc function interpolation
ztrc=sinci(trc,t,tint);
