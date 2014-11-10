function gout=blint(g,t,flow,fhigh,delflow,delfhigh)
 
% gout=blint(g,t,flow,fhigh,delflow,delfhigh)
%
% BLINT performs band limited integration in the frequency domain.
% This is done by spectral division by the first power of frequency
% and an appropriate phase shift. Outside the indicated frequency
% band, the integration filter rolls off as a rapid, smooth gaussian.
%
% g ... digital time series to be integrated
% t ... time coordinate vector for g
% flow ... lowest frequency in passband of integration filter
% fhigh ... highest frequency in passband of integration filter
% delflow ... gaussian rolloff width on low end
%         ******* default min([5 flow/5]) *******
% delfhigh ... gaussian rolloff width on high end
%         ******* default min([(fnyquist-fhigh)/10 10]) *******
%
% G.F. Margrave May 1995
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
%transform (with power of 2 pad)
gp=padpow2(g);
tp=xcoord(t(1),t(2)-t(1),gp);
[G,f]=fftrl(gp,tp);
if(nargin < 6)
	delfhigh= min([f(length(f))-fhigh 10]);
end
if(nargin < 5)
	delflow= min([flow/5 5]);
end
% design integrator
S= ones(size(f))+i*zeros(size(f));
%passband
ind=between(flow,fhigh,f,2);
S(ind) = S(ind)*(-i)./(2*pi*f(ind));
%low end
f1=f(ind(1));
ilow=find(f<f1);
S(ilow)=S(ind(1))*ones(size(ilow));
gw=gauss(f,f1,delflow);
S(ilow)= gw(ilow).*S(ilow);
%high end
f2=f(ind(length(ind)));
ihigh=find(f>f2);
S(ihigh)=S(ind(length(ind)))*ones(size(ihigh));
gw=gauss(f,f2,delfhigh);
S(ihigh)= S(ihigh).*gw(ihigh);
%apply
G=G.*S;
%compute the constant term (assume integration starts at t=0.0)
%c=2*sum(imag(G));
%transform back and truncate
%gout= ifftrl(G,f) + c;
gout= ifftrl(G,f);
gout=gout(1:length(g));
