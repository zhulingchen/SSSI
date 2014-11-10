function [phi,obj]=constphase2(s1,s2,flag)
% CONSTPHASE: find the best constant phase rotation to match one signal to another
%
% [phi,obj]=constphase2(s1,s2,flag)
%
% Given two signals s1 and s2, this function finds the constant phase
% rotation which, when applied to s1, minimizes the L2 norm of the
% difference: s2-phase_rotated(s1) . ("Constant phase" means the phase
% is the same for all frequencies.) 
%
% s1 ... input time series to be rotated
% s2 ... input time series to match to rotated s1
% flag ... 1 means impose the bandwidth of s1 on s2 before determining
%       rotation (see bandwidth_xfer)
%          0 means don't do that
% ************* default = 1 ***********
% phi ... best constant phase angle in degrees. To rotate s1 to look like
%           s2, use s1r=phsrot(s1,phi)
% obj ... the L2 norm of the difference for integer phase rotations from
%           -180 to 179. Phi comes from the minimum of this function.
% NOTE: You can plot obj with: plot(-180:179,obj);xlabel('Angle in degrees')
%
% 
% by G.F. Margrave, Nov. 2013
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

if(nargin<3)
    flag=1;
end

angs=-180:179;
obj=zeros(size(angs));
if(flag==1)
    % impose the bandwidth of the first signal on the second
    s2p=bandwidth_xfer(s1,s2);
else
    s2p=s2;
end
%
for k=1:length(angs)
    s1r=phsrot(s1,angs(k));
    %obj(k)=norm(s2p-s1r);
    obj(k)=sum((s2p-s1r).^2);
end

[om,kmin]=min(obj);
phi=angs(kmin);