function [stp,tau,p]=tptran(seis,t,x,pmin,pmax,dp)
%TPTRAN ... forward tau-p transform (linear trajectories)
%
% [stp,tau,p]=tptran(seis,t,x,pmin,pmax,dp)
%
% The forward tau-p transform is computed in the space-frequency domain.
%
% seis ... input seismic matrix or shot record
% t ... time (row) coordinate for seis
% x ... space (column) coordinate for seis
% NOTE: size(seis,1) must equal size(t,1) and size(seis,2) must equal size(x,2)
% pmin ... minimum slowness value to steer over
% pmax ... maximum slowness value to steer over
% dp ... slowness increment
% ******** default = 0.5*(pmax-pmin)/(length(x)-1) ********
% NOTE: The default for dp gives twice as many p traces as x traces.  The
% corresponding default in itptran gives half as many x traces as p traces.
%
% stp ... tau-p transform of seis
% tau ... row coordinate of stp
% p   ... column coordinate of stp
% 
%
% by G.F. Margrave, September 2014
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

if(size(seis,1)~=size(t,1))
    error('t vector incompatible with seis')
end
if(size(seis,2)~=size(x,2))
    error('x vector incompatible with seis')
end

if(nargin<6)
    dp=0.5*(pmax-pmin)/(length(x)-1);
end

[nt,nx]=size(seis);
nt2=2^nextpow2(nt);
if(nt2>nt)
    seis=[seis;zeros(nt2-nt,nx)];
    tau=((t(2)-t(1))*(0:nt2-1))';
else
    tau=t;
end
%Transform over time
[seisf,f]=fftrl(seis,tau);%seisf is in (x,f) space
p=pmin:dp:pmax;%slowness values to steer over
np=length(p);
stp=zeros(nt2,np);%preallocate space for tau-p transform
%loop over p values
for k=1:length(p)
    dtx=p(k)*x;
    %calculate phase shift operator to flatten linear events
    shiftr=exp(1i*2.*pi*f*dtx);
    %apply the phase shift and then sum (stack) traces
    trcf=sum(seisf.*shiftr,2);
    trcf(end)=real(trcf(end));%make sure Nyquist is real
    stp(:,k)=ifftrl(trcf,f);%inverse transform
end