function [Q,Qrand]=fakeq(vp,rho,Qmin,Qmax,delQpct,vp0,vp1,rho0,rho1,nv,nr)
% FAKEQ: invents a Q log given sonic and density logs
%
% [Q,Qrand]=fakeq(vp,rho,Qmin,Qmax,delQpct,vp0,vp1,rho0,rho1,nv,nr)
% 
%Idea: map Qmin to vp0 and Qmax to vp1 via:
%   Qv=Qmin*((vp-vp0)/(vp1-vp0)).^nv+Qmax*((vp-vp1)/(vp0-vp1)).^nv
% and map Qmin to rho0 and qmax to rho1 via
%   Qr=Qmin*((rho-rho0)/(rho1-rho0)).^nr+Qmax*((rho-rho1)/(rho0-rho1)).^nr
% Then the final Q is given by
%   1./Q = 1./Qv + 1./Qr = (Qv + Qr)./(Qv.*Qr);
%
% WARNING: This is entirely an empirical hunch. No claim of physical
% accuracy is made or implied.
%
% vp ... input velocity 
% rho ... input density
% ***** rho and vp must be vectors of the same size ******
% Qmin ... minimum allowed Q
% Qmax ... maximum allowed Q
% delQpct ... standard deviation of random Q fluctuations, expressed as a
%       percentage (a number between 0 and 100)
% vp0 ... low-end reference velocity
% vp1 ... high end reference velocity
% rho0 ... low-end reference density
% rho1 ... high-end reference density
% nv ... velocity exponent
% nr ... density exponent
%
% Q ... estimated Q vector (deterministic)
% Qrand ... estimated Q with random fluctuations
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

if(size(vp)~=size(rho))
    error('vp and rho must be the same size');
end

Qv=Qmin*((vp-vp1)/(vp0-vp1)).^nv+Qmax*((vp-vp0)/(vp1-vp0)).^nv;
%search for negative Q's
ind=find(Qv<Qmin);
if(~isempty(ind))
    Qv(ind)=Qmin;
end

Qr=Qmin*((rho-rho1)/(rho0-rho1)).^nr+Qmax*((rho-rho0)/(rho1-rho0)).^nr;
%search for negative Q's
ind=find(Qr<Qmin);
if(~isempty(ind))
    Qr(ind)=Qmin;
end

Q= Qv.*Qr./(Qv+Qr);



%include random fluctuations
Qfluct=randn(size(Q)).*Q*delQpct/100;

Qrand=Q+Qfluct;