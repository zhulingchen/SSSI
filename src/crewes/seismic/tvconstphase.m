function [phs,tout]=tvconstphase(s1,s2,t,twin,tinc,flag)
% TVCONSTPHASE: estimates apparent temporally local constant phase rotations 
%
% [phs,tout]=tvconstphase(s1,s2,t,twin,tinc,flag)
% 
% The seismic trace s1 is localized in time with a Gaussian window and then
% the constant phase rotation which best matches the localized trace to s2
% (with the same window applied) is computed (see constphase). This process
% is repeated until all specified times are analyzed.
%
% s1= input trace to be analyzed
% s2= reference trace. Constant phase rotations are w.r.t. this trace
% t= time coordinate vector for s1
% twin= width (seconds) of the Gaussian window
% tinc= temporal shift (seconds) between windows
% flag ... 1 means impose the bandwidth of s1 on s2 before determining
%       rotation (done independently for each window) (see bandwidth_xfer)
%          0 means don't do that
% ************* default = 1 ***********
% phs= apparent constant phase rotations in each window
% tout= window center times. Same size as phs
%
% by G.F. Margrave, Sept 2005
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

if(nargin<6)
    flag=1;
end

tmin=t(1);
t=t-tmin;
% dt=t(2)-t(1);
%determine number of windows. tinc will be adjusted to make the
% last window precisely centered on tmax
tmax=t(end);
nwin=tmax/tinc+1; %this will generally be fractional
nwin=round(nwin);
tinc=tmax/(nwin-1); %redefine tinc
tout=zeros(nwin,1);
phs=tout;
for k=1:nwin
    %build the gaussian
    tnot=(k-1)*tinc;
    tout(k)=tnot;
    gwin=exp(-((t-tnot)/twin).^2)/(sqrt(pi)*twin/tinc);
    %window and measure phase
    s1w=s1.*gwin;
    s2w=s2.*gwin;
    phs(k)=constphase2(s1w,s2w,flag);
end
tout=tout+tmin;