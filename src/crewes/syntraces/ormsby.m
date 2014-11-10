function [w,t]=ormsby(f1,f2,f3,f4,tlen,dt)
% ORMSBY: Creates an Ormsby bandpass filter
%
% [w,t]=ormsby(f1,f2,f3,f4,tlen,dt)
%
% Make a ormsby wavelet given the four frequency parameters
% f1 = low frequency stop
% f2 = lowest frequency to pass unattenuated
% f3 = highest frequency to pass attenuated`
% f4 = high frequency stop
% tlen = wavelet length in seconds. Wavelet will go from
%	-tlen/2 to tlen/2
% dt = wavelet sample rate in seconds
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

% make the time coordinate vector
    nt=tlen/dt+1;
    nt2=floor(nt/2);
    t=dt*(-nt2:nt-nt2-1)'; %ensures a sample at t=0

%make the wavelet

	w= (pi*f4^2) * (sinque(pi*f4*t)).^2/(f4-f3);
	w= w- (pi*f3^2) * (sinque(pi*f3*t)).^2/(f4-f3);
	w= w- (pi*f2^2) * (sinque(pi*f2*t)).^2/(f2-f1);
	w= w+ (pi*f1^2) * (sinque(pi*f1*t)).^2/(f2-f1);
 
%normalize
    w=wavenorm(w,t,2);
