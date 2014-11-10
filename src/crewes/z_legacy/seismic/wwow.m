function [wavelet,tw,pseudo]=wwow(trin,t,tleng)
% [wavelet,tw]=wwow(trin,t,tleng)
%
% WWOW performs the 'wwow' method of wavelet extraction which
% is identical to a single channel BOOT. The final wavelet is 
% determined via a time domain match filter.
%
% trin= input trace
% t= time coordinate for input trace
% tleng= maximum length of output wavelet
% wavelet= output wavelet
% tw= time coordinate for wavelet
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
 %trin=trin(:)';
% get the hilbert mask
 [mask,htrin]=hmask(trin);
% find the constant phase rotation
  hsum=maxima(htrin,mask);
 top= 2.*sum(real(hsum).*imag(hsum));
 bot= sum(imag(hsum).^2-real(hsum).^2);
 theta=.5*atan2(top,bot);
% test theta and theta+pi/2
 trot1=phsrot(trin,theta*180/pi)';
 trot2=phsrot(trin,(theta+pi/2)*180/pi)';
 rms1=norm(trot1.*mask);
 rms2=norm(trot2.*mask);
 if rms2>rms1, theta=theta+pi/2; trot1=trot2; end
% compute the pseudo reflectivity
 pseudo=mask.*trot1;
 wavelet=theta*180/pi;
%
% match filter for the wavelet
%
  dt=t(2)-t(1);
  wavelet=match(pseudo,trin,t,tleng,0);
  tw=xcoord(0.,dt,wavelet);
 
