function wavemin=tomin(wavelet,stab)
% wavemin= tomin(wavelet,stab)
% wavemin= tomin(wavelet)
%
% TOMIN uses the Levinson algorithm to convert an arbitrary
% waveform to its minimum phase equivalent.
%
% wavelet= input waveform to be converted
% stab= stab factor to be used expressed as a fraction of the
%       zero lag of the aurocorellogram
%   *********** default = .0001 ************
% wavemin= output minimum phase equivalent of wavelet
%
% The wavelet is autocorrelated and input to LEVREC.
% The output minimum phase wavelet is the frequency domain 
% inverse of the output from LEVREC.
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
if(nargin<2)
   stab=.0001;
 end
% generate the auto
  nlags=length(wavelet);
  a=auto(wavelet,nlags,0);
  a(1)=a(1)*(1.+stab);
  a= a/a(1);
% run this through Levinson
  b=[1.0 zeros(1,nlags-1)]';
  winv=levrec(a(1:nlags),b);
% invert the wavelet
  wavemin=real(ifft(1. ./(fft(winv'))));
% now normalize the wavelet
  powrms=sqrt(wavemin*wavemin');
  wavemin=wavemin/powrms;
 
 		
    
