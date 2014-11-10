function [freqloc,freqins]=fomelfreq(s,t,tsmo,lambda)
% FOMELFREQ ... calulate a local frequency based on a method by S. Fomel
%
% [freqloc,freqins]=fomelfreq(s,t,tsmo,lambda)
%
% Method: See "Local seismic attributes", by S. Fomel, Geophysics, 2007
%
% s = signal
% t = time coordinate for s
% tsmo = halfwidth of triangle smoother (seconds)
% lambda = stabilization constant
%
% freqloc ... Fomel's local frequency
% freqins ... classic instantaneous frequency
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

s=s(:);
dt=t(2)-t(1);
%differentiate signal by central finite difference
sprime=zeros(size(s));
sprime(2:end-1)=(s(3:end)-s(1:end-2))/(4*dt*pi);
%hilbert transforms
sh=imag(hilbert(s));
sprimeh2=imag(hilbert(sprime));
sprimeh=zeros(size(s));
sprimeh(2:end-1)=(sh(3:end)-sh(1:end-2))/(4*dt*pi);
%numerator of equation 2
top=s.*sprimeh-sprime.*sh;
%denominator of equation 3
bot=s.^2+sh.^2;
%big D
D=diag(bot);
%make the triangle smoothing matrix
nsmo=round(tsmo/dt);
tri=triang(2*nsmo+1);
S=convmtx(tri/sum(tri),length(s));
S(1:nsmo,:)=[];
S(end-nsmo+1:end,:)=[];
%now build the denominator of equation 7
Bot=lambda^2*eye(length(s))+S*(D-lambda^2*eye(length(s)));
%calculate both frequencies
freqloc=Bot\S*top;
freqins=top./bot;



