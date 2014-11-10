function [quimp,tq,b]=einar(Q,x,cnot,dt,tmax)
% [quimp,tq]=einar(Q,x,cnot,dt,tmax)
%
% Generate a constant Q impulse response
%
% Q ... Q
% x ... distance traveled
% cnot ... high frequency phase velocity
% dt ... time sample rate
% tmax ... maximum time output
% quimp ... impulse response (in retarded time)
% tq ... time vector for quimp
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
%fnot= f(length(f));
tq=0:dt:tmax;
%Find a power of 2 thats at least 4 times the length of tq
nf=length(tq);
while(nf < 4*length(tq) )
	nf = 2^nextpow2(nf+1);
end
%build frequency vector
fnyq=1/(2*dt);
df=1/((nf-1)*dt);
f=0:df:2*fnyq;
ind=nf/2 + 1;
f=[f(1:ind) -f(ind-1:-1:2)];
%amplitude spectrum
B= exp(-x*pi*abs(f)/(Q*cnot));
%phase velocity
cw=cnot*(1+log(abs(f(2:length(f)))/fnyq)/(pi*Q));
cw=[cw(1) cw];
%complex spectrum
B = B.*exp(i*2*pi*f.*(x/cnot - x./cw));
%to time domain
b=real(ifft(B));
%apply window
quimp=b(1:length(tq)).*mwhalf(length(tq));
