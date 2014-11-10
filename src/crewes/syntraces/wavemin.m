function [w,tw]=wavemin(dt,fdom,tlength,m,stab)
% WAVEMIN: creates a minimum phase wavelet for impulsive sources
%
% [w,tw]=wavemin(dt,fdom,tlength,m,stab)
% [w,tw]=wavemin(dt,fdom,tlength,m)
% [w,tw]=wavemin(dt,fdom,tlength)
% [w,tw]=wavemin(dt,fdom)
% [w,tw]=wavemin(dt)
% 
% WAVEMIN returns a minimum phase wavelet which simulates a possible noise
% free impulsive source. Function uses TNTAMP to generate a plausible power
% spectrum. This is then inverse Fourier transformed into an
% autocorrelation which is fed into the Levinson recursion (LEVREC).
% Finally, the result from the Levinson recursion is inverted (by frequency
% domain division) to get the final wavelet.  The wavelet is normalized
% such that a sine wave at the dominant frequency is passed with unit
% amplification.
%
% dt= desired temporal sample rate (seconds)
% fdom= dominant frequency in Hz
%  ******** default: 15 Hz *******
% tlength= w length in seconds 
%  ******** default: 127*dt (ie a power of 2) *******
% m = exponent controlling spectral shape. See tntamp for a description.
%       Larger values give more rapid falloff of high frequencies. The
%       default is a moderate decay. A value like 10 or 12 will be extreme.
%       Put another way, a low value for m (like 2) gives a wavelet with a
%       very broad band that is easy to deconvolve. A high number (like 7)
%       is very narrow band and may not even be truly minimum phase. See
%       the example below.
% ************ default 4 ************
% stab = white noise factor. In order to stabilize the Levinson recursion,
%       the zero-lag autocorrelation is multiplied by 1+stab. This matters
%       for values of m larger than 2.
% ************ default = .000001 **********
%
% Example: test different m values
% Copy/paste this to the command line to see what m does.
%
% [w2,tw]=wavemin(.001,30,.2,3);
% [w3,tw]=wavemin(.001,30,.2,2);
% [w4,tw]=wavemin(.001,30,.2,4);
% [w5,tw]=wavemin(.001,30,.2,5);
% [w6,tw]=wavemin(.001,30,.2,6);
% [w7,tw]=wavemin(.001,30,.2,7);
% figure
% dbspec(tw,[w2 w3 w4 w5 w6 w7])
% legend('m=2','m=3','m=4','m=5','m=6','m=7')
% 
% by G.F. Margrave, May 1991, 2013
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

if(nargin<5); stab=.000001; end 
if(nargin<4); m=4; end
 if(nargin<3)
   tlength=127.*dt;
 end
 if(nargin<2)
   fdom=15.; 
 end
% create a time vector
  nt=round(2.*tlength/dt)+1;
  nt=2.^nextpow2(nt);
  tmax=dt*(nt-1);
  tw= 0.:dt:tmax;
% create a frequency vector
  fnyq=1./(2*(tw(2)-tw(1)));
  f=linspace(0.,fnyq,length(tw)/2+1);
% create the power spectrum
  tmp=tntamp(fdom,f,m);
  %tmp=(tmp/max(tmp)).^4;
  powspec=tmp.^2;
% create the autocorrelation
  auto=ifftrl(powspec,f);
  auto(1)=auto(1)*(1+stab);
% run this through Levinson
  nlags=round(tlength/dt)+1;
  b=[1.0 zeros(1,nlags-1)]';
  winv=levrec(auto(1:nlags),b);
% invert the winv
  w=real(ifft(1. ./(fft(winv))));
  tw=(dt*(0:length(w)-1))'; 
% now normalize the w
  w=wavenorm(w,tw,2);
