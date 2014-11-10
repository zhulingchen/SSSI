function winv=toinvf(wavelet,stab,n)
% TOINVF FFT based waveform inversion with low-pass zero-phase filter
% winv=toinvf(wavelet,stab,n)
% winv=toinvf(wavelet,stab)
% winv= toinvf(wavelet)
%
% TOINVF uses FFT to invert an arbitrary waveform over its
% entire frequency band except for a zero phase high frequency
% attenuation. A stab factor is supplied to avoid
% instabilities caused by division by small numbers.
%
% wavelet= input waveform to be converted
% stab= stab factor to be used expressed as a fraction of the
%       peak of the amplitude spectrum 
%   *********** default = .001 ************
% n= length of inverse (number of points)
%  ********** default= length(wavelet) ***********
% winv= output inverse wavelet. Winv will be the same length as
%       wmin. It is occaisionally advisable to pad wmin with zeros
%       before inverting.
%
% note: this inverse is not usually causal, use CONVZ to apply
% it (default time zero to length(winv)/2)
%
% by G.F. Margrave, June 1991
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

 if nargin<3, n=length(wavelet); end
 if(nargin<2)
   stab=.001;
 end
% prepare the wavelets
 lw=max([length(wavelet),n]);
 wavelet=padpow2(wavelet,1);
 while (length(wavelet)/lw)<1.5
  wavelet=padpow2(wavelet,1);
 end
  imp=impulse(wavelet);
% invert the wavelet
  Wlet=fft(wavelet);
  Imp=fft(imp);
  bground=zeros(size(Wlet))+stab*max(abs(Wlet));
  indicies= find(abs(Wlet)<bground);
  Wlet(indicies)=bground(indicies);
% invert
  Winv=Imp./Wlet;
% band limit it
  mw=mwindow(length(Winv)+1,10);
  mw=fftshift(mw(1:length(Winv)));
  Winv=Winv.*mw;
% back to time
  winv=real(ifft(Winv));
  n2=length(winv)/2;
  winv1=winv(round(n2-n/2+1):round(n2+n/2)).*hanning(n)';
  winv=winv1;

 




 		
    
