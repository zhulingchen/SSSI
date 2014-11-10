function [m,a,rhonew]=topsgardner(vel,rho,t,tops)
% TOPSGARDENER defines a set of Gardner parameters that can be used to 
%     calculate a density log from a velocity log.  This algoritim 
%     defines a new set of paramters for each lithology interval as defined
%     in tops. This is meant to be applied in time but can be applied in 
%     depth as long as all input parameters are in depth. 
%
%     The output parameters can be applied using the following code:
%             rhonew=m.*vel.^a;
%
% [m,a,rhonew]=topsgardner(vel,rho,t,tops)
%
% vel    = seismic data matrix
% rho    = time vector
% t      = well log to match amplitudes to
% tops   = time corresponding to the end of the log
%
% m      = velocity - density scaling factor
% a      = velocity - density exponent factor
% rhonew = the new density vector as calculated using the topsgardner
%            parameters
%
% by  H.J.E. Lloyd December 2012
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
vel=vel(:);
rho=rho(:);
t=t(:);
tops=sort(tops(:));
if tops(1)~=0
    tops=[0;tops];
end

if tops(end)~=t(end)
    tops=[tops;t(end)];
end
tops=sort(tops(:));
tmin=t(1);
t=t-tmin;
%determine number of windows. tinc will be adjusted to make the
% last window precisely centered on tmax
a=zeros(size(t));
m=a;
for k=1:length(tops)-1
    ind=near(t,tops(k)):near(t,tops(k+1));
pp=polyfit(log(vel(ind)),log(rho(ind)),1);
aout=pp(1);
mout=mean(rho(ind)./(vel(ind).^aout));
m(ind)=mout;
a(ind)=aout;
end
rhonew=m.*vel.^a;




