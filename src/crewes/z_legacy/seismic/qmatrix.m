function qmat=qmatrix(q,t,w,tw)
% qmat=qmatrix(q,t,w,tw)
%
% Generates a matrix which applies a forward q filter. The qmatrix is
% also bandlimited by a stationary "source signature" waveform.
%
% q ... value of Q
% t ... vector of times to generate q response at
% w ... stationary source waveform
% tw ... time vector for w
% qmat ... length(t) by length(t) matrix which applies a forward Q
%		filter.
% 		If r is a column vector of reflection coefficients, then
%		s=qmat*r; is a simple synthetic seismogram with attenuation.
%
% G.F. Margrave, July 1996, The CREWES Project
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
v=2000;
x=v*t;
nt=length(t);
tmax=max(t);
dt=t(2)-t(1);
qmat=zeros(nt,nt);
if(nargin>2)
	izero=near(tw,0);
end
for k=1:nt
	tmp=einar(q,x(k),v,dt,tmax);
	nuse=round((tmax-t(k))/dt) +1;
	
	%include stationary wavelet
	if(nargin>2)
		tmp=convz(tmp,w,izero,nt,0);
	end
	
	qmat(k:nt,k)=tmp(1:nuse)';
end
%wconv=convmtx(w(:),nt);
%qmat=wconv(1:nt,:)*qmat;
	
	
