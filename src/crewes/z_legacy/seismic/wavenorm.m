function w=wavenorm(win)
% w=wavenorm(win)
%
% win= input wavelet
%
% WAVENORM normalizes a wavelet in the sense that the a sinusoidal
% signal at the peak of the wavelets passband will be unchanged in
% amplitude if the wavelet is convolved with said signal.
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
%determine the peak frequency
%fake a time axis
t=xcoord(0,.002,length(win));
%compute spectrum
[s,f]=fftrl(win,t);
s=abs(s);
%find the max of s
maxs=max(s);
ifreq=find(s==maxs);
fnot=f(ifreq);
%make sure there is only 1 fnot
fnot=fnot(length(fnot));
%generate a reference sinusoid
ref=sin(2*pi*fnot*t);
%convolve
c=conv(ref,win);
%find the max of c
cmax=max(c);
%normalize
w=win/cmax;
