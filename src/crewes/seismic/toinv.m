function winv=toinv(wavelet,stab,n)
% TOINV inverts an arbitrary waveform using match filtering
% winv= toinv(wavelet,stab,n)
% winv= toinv(wavelet,stab)
% winv=toinv(wavelet)
%
% TOINV uses match filtering to invert an arbitrary waveform.
% For a band limited inverse, simply follow this by a filter. 
% A stab factor is supplied to avoid instabilities caused by 
% division by small numbers.
%
% Note: This inverse is always causal and should be applied
% with CONV or CONVM. Sometimes a better, though 
% non-causal, inverse is obtained by TOINVF.
%
% wavelet= input waveform to be converted
% stab= stab factor to be used expressed as a fraction of the
%       peak of the amplitude spectrum 
%       **** default = .001 ****
% n= length of inverse operator
%    **** default = length(wavelet) ****
% winv= output inverse wavelet
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

 if nargin<3, n=length(wavelet);end
 if(nargin<2)
   stab=.001;
 end
 Wlet=fft(wavelet);
 bground=stab*max(abs(Wlet));
 wavelet=real(ifft(Wlet+bground));
 winv=zeros(1,round(n));
 tw=xcoord(0.,.004,winv);
 winv=match(wavelet,impulse(wavelet,1),tw,max(tw),1);

 




 		
    
