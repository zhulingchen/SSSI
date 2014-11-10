function [seis,t,x]=itptran(stp,tau,p,xmin,xmax,dx)
%ITPTRAN ... inverse tau-p transform (linear trajectories)
%
% [seis,t,x]=itptran(stp,tau,p,xmin,xmax,dx)
%
% The inverse tau-p transform is computed in the space-frequency domain.
% The algorithm is filtered back-projection.
%
% stp ... input tau-p transform
% tau ... time (row) coordinate for stp
% p ... slowness (column) coordinate for stp
% NOTE: size(stp,1) must equal size(tau,1) and size(stp,2) must equal size(p,2)
% xmin ... minimum x value for output
% pmax ... maximum x value for output
% dx ... spatial increment
% ******** default = 2*(xmax-xmin)/(length(p)-1) ********
% NOTE: The default for dx gives have as many x traces as p traces.  The
% corresponding default in tptran gives twice as many p traces as x traces.
%
% seis ... output seismic gather
% t ... row coordinate of seis
% x   ... column (space) coordinate of seis
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

if(size(stp,1)~=size(tau,1))
    error('tau vector incompatible with stp')
end
if(size(stp,2)~=size(p,2))
    error('p vector incompatible with stp')
end

if(nargin<6)
    dx=2*(xmax-xmin)/(length(p)-1);
end

[nt,np]=size(stp);
nt2=2^nextpow2(nt);
if(nt2>nt)
    stp=[stp;zeros(nt2-nt,np)];
    t=(tau(2)-tau(1))*(0:nt2-1);
else
    t=tau;
end
[stpf,f]=fftrl(stp,tau);
x=xmin:dx:xmax;
nx=length(x);
seis=zeros(nt2,nx);
%loop over x values
for k=1:length(x)
    dtx=x(k)*p;
    %calculate phase shift operator to flatten linear events
    shiftr=exp(-1i*2.*pi*f*dtx);
    %apply the phase shift and then sum (stack) traces
    trcf=f.*sum(stpf.*shiftr,2);
    trcf(end)=real(trcf(end));%make sure Nyquist is real
    seis(:,k)=ifftrl(trcf,f);%inverse transform
end