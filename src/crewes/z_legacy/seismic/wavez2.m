function [wavelet,twave]=wavez2(f1,f2,tlength,t)
% [wavelet,twave] = wavez2(f1,f2,tlength,t)
%
% returns a zero phase wavelet of passband as specified
%
% f1= low end of frequency bandpass (Hz)
% f2= high end of frequency bandpass (Hz)
% tlength=temporal length of wavelet
% t= time sample vector with at least two entries 
%     (provides time sample rate)
% wavelet = returned zero phase wavelet
% twave= returned time coordinate vector for wavelet
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
 n=round(tlength/(t(2)-t(1)))+1;
 fnyq= .5/(t(2)-t(1));
 wavelet= fir1(n-1,[f1/fnyq,f2/fnyq]);
 twave=linspace(-tlength/2,tlength/2,n);
% normalize the wavelet
% generate a refenence sinusoid at the dominant frequency
 refwave=sin(2*pi*fdom*twave);
 reftest=convm(refwave,wavelet);
 fact=max(refwave)/max(reftest);
 wavelet=wavelet*fact;
