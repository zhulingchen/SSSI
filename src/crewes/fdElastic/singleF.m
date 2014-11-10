function [wave,tw] = singleF(Dt,freq,cycles)
%[wave,tw] = singleF(Dt,freq,cycles)
%Design (near) single frequency zero phase sweep (tapered)
% Dt     = The time sample rate in seconds
% freq   = The design frequency
% cycles = The basic number of cycles on each side of the centre
% wave   = The wavelet amplitudes in sequence
% tw     = The wavelet time series in sequence
%
% P.M. Manning, Dec 2011
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
%Dt = .0005; %.002;          %Depends on frequency
%freq = 30; %60; %30;      %Frequency in Hz
disp(cycles)
top = 1.25*cycles/freq;
nHalf = ceil(top/Dt);
%nFilt = nHalf*2-1;
high = (nHalf-1)*Dt;
tw = (-high:Dt:high);
signal = cos(tw*freq*2*pi);
A = 1/sqrt(2*pi);
%fact = freq*3.4/(1.25*cycles);
fact = freq*4.5/(1.25*cycles);      %Tapers wavelet end
power = ((fact*tw).^2)*0.5;
env = A*exp(-power);
wave = signal.*env;
%plot(tw,signal,tw,wave)
%plot(tw,env)
disp(wave(1:5))
