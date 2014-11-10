function [wzero,tw]=tozero(wavelet,twin)

% [wzero,tw]=tozero(wavelet,twin)
% wzero=tozero(wavelet)
%
% TOZERO uses FFT to capture the amplitude spectrum of an 
% arbitrary and outputs an 'all pass' equivalent.
%
% wavelet= input waveform to be converted
% twin= input time coordinate
% wzero= output all pass wavelet
% tw= output time coordinate. twin must be supplied for this
%     to be computed
%
% note: this result is never causal, use CONVZ to apply
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

if nargin<2, twin=xcoord(0.,.002,wavelet); end
 lw=length(wavelet);
 [rw,cw]=size(wavelet);
 wavelet=wavelet(:);
 wavelet=padpow2(wavelet);
% fft
 Wlet=fft(wavelet);
 Imp=fft(impulse(wavelet));
 A=abs(Wlet);
 wzero=real(ifft(A.*Imp));
% grab central lw samples
 lw2=length(wavelet);
 wzero=wzero( round(lw2/2-lw/2+1):round(lw2/2+lw/2) ).*mwindow(lw);
% make time vector
 if nargin==2
   n=floor(length(wzero)/2);
   dt=twin(2)-twin(1);
   tw=xcoord(-(n-1)*dt,dt,wzero);
 end
 if(rw==1)wzero=wzero';end


 




 		
    
