function [w,tw,wd]=wavez(dt,fdom,tlength,m)
% WAVEZ: creates a zero phase w with a dominant frequency
%
% [w,tw]=wavez(dt,fdom,tlength,m)
% [w,tw]=wavez(dt,fdom,tlength)
% [w,tw]=wavez(dt,fdom) 
% [w,tw]=wavez(dt) 
% 
% WAVEZ returns a zero phase w with a realistic amplitude spectrum
%
% dt= desired temporal sample rate
% fdom= dominant frequency in Hz (default: 15 Hz)
% tlength= w length in seconds (default: 127*dt 
%                                     (ie a power of 2))
% m = exponent controlling spectral shape. See tntamp for a description
% ************ default 2 ************
% 
% 
% by G.F. Margrave, July 1991
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

if(nargin<4) m=2; end
 if(nargin<3)
   tlength=127.*dt;
 end
 if(nargin<2)
   fdom=15.; 
 end

% create a time vector
  nt=tlength/dt+1;
  nt=2.^nextpow2(nt);
  tmax=dt*(nt)/2;
  tw=-tmax+dt*(0:nt-1)';;
  fnyq=1./(2*(tw(2)-tw(1)));
  f=linspace(0.,fnyq,length(tw)/2+1)';
% create the amplitude spectrum spectrum
  aspec=tntamp(fdom,f,m);
% inverse transform
  w=ifftrl(aspec,f);
  w=fftshift(w);
  wd=ifftrl(i*f.*aspec,f);
  wd=fftshift(wd);
% make sure we don't return the Fourier pad if it was not asked for
  ind=between(-tlength/2,tlength/2,tw,2);
  w=w(ind);
  tw=tw(ind);
  wd=wd(ind);
% now normalize the w
  w=wavenorm(w,tw,2);