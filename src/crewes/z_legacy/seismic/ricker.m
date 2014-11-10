function [wavelet,tw]=ricker(dt,fdom,tlength)
% [wavelet,tw]=ricker(dt,fdom,tlength)
% [wavelet,tw]=ricker(dt,fdom) 
% [wavelet,tw]=ricker(dt) 
% 
% RICKER returns a Ricker wavelet.
%
% dt= desired temporal sample rate
% fdom= dominant frequency in Hz (default: 15 Hz)
% tlength= wavelet length in seconds (default: 127*dt 
%                                     (ie a power of 2))
% 
% The wavelet is generated from an analog expression for a 
% Ricker wavelet.
% 
% by G.F. Margrave, May 1991
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
if(nargin<3)
   tlength=127.*dt;
 end
 if(nargin<2)
   fdom=15.; 
 end
% create a time vector
  nt=round(tlength/dt)+1;
  tmin=-dt*round(nt/2);
  tmax=-tmin-dt;
  tw= tmin:dt:tmax;
% create the wavelet
  pf=pi^2*fdom^2;
  wavelet=(1-2.*pf*tw.^2).*exp(-pf*tw.^2);
% normalize
% generate a refenence sinusoid at the dominant frequency
  refwave=sin(2*pi*fdom*tw);
  reftest=convz(refwave,wavelet);
  fact=max(refwave)/max(reftest);
	wavelet=wavelet*fact;
