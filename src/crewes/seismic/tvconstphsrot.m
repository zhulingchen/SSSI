function [trout,phi,tphi]=tvconstphsrot(trin,tref,t,twin,tinc)
% TVCONSTPHSROT Constant-phase rotate a trace for gaussian time windows
% [trout,phi,tphi]=tvconstphsrot(trin,tref,t,twin,tinc)
%
% TVCONSTPHSROT performs a time varient constant phase rotation of the input trace
% through an vector of angles (phi) in degrees.
%
% trin= input trace
% tref= reference trace. Constant phase rotations are w.r.t. this trace
% t= time coordinate vector for trin
% twin= width (seconds) of the Gaussian window
% tinc= temporal shift (seconds) between windows
%
% trout= phase rotated trace
% phi=   Phase rotation angles
% tphi=  time vector for phi

% by G.F. Margrave and H.J.E. Lloyd Jan 2012
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
tmin=t(1);
t=t-tmin;
dt=t(2)-t(1);
%determine number of windows. tinc will be adjusted to make the
% last window precisely centered on tmax
tmax=t(end);
nwin=tmax/tinc+1; %this will generally be fractional
nwin=round(nwin);
tinc=tmax/(nwin-1); %redefine tinc
tout=zeros(nwin,1);
tout=(0:nwin-1)'.*tinc(:)+tmin;
[gwin,norm_factor,tinc2,nwin]=gaussian_upou(t,tmin,twin,tinc);
phi=zeros(size(tout));
tphi=tout;
trmat=zeros(length(trin),nwin);
%loop over windows
wbar=waitbar(0,'Please Wait...');
kk=0;
for k=1:nwin
    %build the gaussian
    tnot=(k-1)*tinc;
    if tnot<=max(t)
    gwin=gaussian_upou(t,tnot,twin,tinc,norm_factor);
    %build the gaussian
    tnot=(k-1)*tinc;
    tout(k)=tnot;
    gwin=exp(-((t-tnot)/twin).^2)/(sqrt(pi)*twin/tinc);
    %window and measure phase
    s1w=trin.*gwin;
    s2w=tref.*gwin;
    phi(k)=constphase(s1w,s2w);
    waitbar(k/nwin,wbar);
    kk=kk+1;
    else
        break
    end
end
phi=phi(1:kk);
tphi=tout(1:kk);
delete(wbar);
phit=ppval(pchip(tphi,phi),t);
trout=trin.*cosd(phit)+phsrot(trin,90).*sind(phit);
