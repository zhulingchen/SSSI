function trout=calfil(trin,t,aec_length,flow,fhigh,nfilt)
% trout=calfil(trin,t,aec_length,flow,fhigh,nfilt)
% 
% CALFIL emulates Chevron's time variant spectral whitening 
% program of the same name. A set of gaussian filters are
% designed to slice a frequency range in such a way that the
% linear sum of slices returns the original spectrum. Each slice
% is then ifft'd into the time domain and aec'd. Finally, the 
% slices are summed in the time domain.
%
% trin= input trace 
% t= time coordinate vector for trin
% aec_length= length (seconds) of the aec operator
% flow= lowest frequency to be whitened (Hz)
% fhigh= highest frequency to be whitened (Hz)
% nfilt= number of Gaussian filter slices
% trout= output trace
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
  
%forward fft
 m=length(trin);
 trout=zeros(size(trin));
 trin=padpow2(trin);
 t=xcoord(t(1),t(2)-t(1),length(trin));
 [spec,f]=fftrl(trin,t);
% loop over slices
 fwidth=(fhigh-flow)/(nfilt-1);
 fnot=flow;
 for n=1:nfilt
   g=gauss(f,flow+(n-1)*fwidth,fwidth); % make gaussian
   Slice=spec.*g; % make slice
   slice=ifftrl(Slice,f); % inverse transform
   slice=aec(slice,t,aec_length); % aec
   trout=trout+slice(1:m); % accumulate slices
 end
 trout=balans(trout,trin);
    
 
  
  
