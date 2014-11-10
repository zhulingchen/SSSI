function aout=hyperbola(x,t,xnot,tnot,v,fdom)
% aout=hyperbola(x,t,xnot,tnot,v,fdom)
%
% Generate a set of unit spikes at an approximate hyperbolic trajectory.
% A ricker wavelet of dominant frequency  fdom is then convolved.
%
% x ... vector of x coordinates of traces
% t ... vector of time coordinates of samples
% xnot ... x location of hyperbola
% tnot ... t location of hyperbola
% v ... velocity
% fdom ... ricker wavelet dominant frequency
% aout ... output matrix of length(t) rows by length(x) columns
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
nsamp=length(t);
ntr=length(x);
%make the wavelet
tlen=.2*(t(nsamp)-t(1));
dt=t(2)-t(1);
[w,tw]=ricker(dt,fdom,tlen);
nw=length(w);
aout=zeros(nsamp+nw,ntr);
v=v/2;
for k=1:ntr
	ispike = round(1+ sqrt( tnot*tnot + ((x(k)-xnot)/v)^2 )/dt);
	if( ispike< (nsamp+nw) )
		aout(ispike,k)=1.0;
		aout(:,k) = convz(aout(:,k),w);
	end
end
aout=aout(1:nsamp,:);
