function [phiout,seisout]=tvlsphsrot(seis,x,t,wells,xwells,twin,tinc,avgflg)
% TVPHSROT Constant-phase rotate a trace
% trout=tvphsrot(trin,phi,t,twin,tinc)
%
% TVPHSROT performs a time varient constant phase rotation of the input trace
% through an vector of angles (phi) in degrees.
%
% trin= input trace
% phi= vector of phase rotation angles in degrees created by tvconstphase
% t= time coordinate vector for trin
% twin= width (seconds) of the Gaussian window
% tinc= temporal shift (seconds) between windows
%
% trout= phase rotated trace

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
phs=tout;
winmat=zeros(length(t),nwin);
for k=1:nwin
    %build the gaussian
    tnot=(k-1)*tinc;
    %build the gaussian
    tnot=(k-1)*tinc;
    tout(k)=tnot;
    gwin=exp(-((t-tnot)/twin).^2)/(sqrt(pi)*twin/tinc);
    winmat(:,k)=gwin(:);
end
trcs=zeros(length(t),length(xwells));
for k=1:size(wells,2)
    ind=near(x,xwells(k));
    ind=ind-avgflg:ind+avgflg;
    if any(ind)<0
        ind=ind(ind>0);
    end
    if any(ind)>length(x)
        ind=ind(ind<length(x));
    end
trcs(:,k)=mean(seis(:,ind),2);
end
tout=(0:nwin-1)*tinc+tmin;
%loop over windows
wbar=waitbar(0,'Please Wait...');
for k=1:nwin
    %build the super trace 
    sw1=[];
    st1=[];
    for kk=1:size(wells,2)
      sw1=[sw1(:);wells(:,kk).*winmat(:,k)];
      st1=[st1(:);trcs(:,kk).*winmat(:,k)];
    end
    phi(k)=constphase(st1,sw1);
    waitbar(k/nwin,wbar);
end
delete(wbar);

phiout=ppval(pchip(tout(:),phi(:)),t(:));



for k=1:size(x,2)
    seisp=phsrot(seis(:,k),90);
seisout(:,k)=seis(:,k).*cosd(phiout(:))+seisp.*sind(phiout(:));
end
end