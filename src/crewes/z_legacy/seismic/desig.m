function trout=desig(trin,wavelet,stab,flag)
%
% trout= desig(trin,wavelet,stab,flag)
% trout= desig(trin,wavelet,stab)
% trout= desig(trin,wavelet)
%
% DESIG uses frequency domain methods to deconvolve a wavelet out
% of trin. A stab factor is added to the spectrum of wavelet which 
% is then inverted over its entire frequency band (except for a 
% zero phase high frequency rolloff). The spectrum of trin is 
% then multiplied by this operator spectrum and inverse transformed.
% A sufficient zero pad is automatically added to trin.
%
% trin= input trace to be de-signatured
% wavelet= input waveform (the known signatur)
% stab= stab factor to be used expressed as a fraction of the
%       peak of the amplitude spectrum.  
%   *********** default = .001 ************
% Alternate: If stab > 1, it is interpreted as being expressed
%   in dbdown. Thus stab=80 is equivalent to stab=.0001
% flag= 0 ... wavelet is not causal
%     = 1 ... wavelet is causal
% ************ default=1 ***********
%
% trout= output desig'd trace.
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
if nargin<4, flag=1; end
 if(nargin<3)
   stab=.001;
 end
if stab>1
   stab=10^(-stab/20);
end
% prepare the vectors
 lt=length(trin);
 lw=length(wavelet);
 trin=padpow2(trin);
 while length(trin)<2(lt+lw)
   trin=padpow2(trin);
 end
 wavelet=pad(wavelet,trin);
% invert the wavelet
  Wlet=fft(wavelet);
  bground=zeros(size(Wlet))+stab*max(abs(Wlet));
  indicies= find(abs(Wlet)<bground);
  Wlet(indicies)=bground(indicies);
  Imp=ones(size(Wlet));
  if flag==0
    Imp=fft(pad(impulse(zeros(1,lw)),wavelet));
  end
% invert
  Winv=Imp./Wlet;
  % band limit it
  mw=mwindow(length(Winv)+1,10);
  mw=fftshift(mw(1:length(Winv)));
  Winv=Winv.*mw;
% apply to trin
  Trin=fft(trin);
  Trout= Trin.*Winv;
% back to time
  trout=real(ifft(Trout));
  trout1=trout(1:lt).*mwindow(lt,5);
   trout=trout1;
 
 		
    
