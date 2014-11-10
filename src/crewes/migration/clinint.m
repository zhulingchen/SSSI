function vint = clinint(f,v,fint)
% CLININT: Complex-valued linear interpolation (used by FKMIG)
%
% vint = clinint(f,v,fint)
%
% CLININT does complex linear interpolation of a function assumed to the
% the t->f transform of a time domain function defined from 0 to some tmax.
% In order to avoid spectral artifacts, the interpolation weights are 
% complex (i.e. have phase terms). These weights may be derived by 
% considering sinc function interpolation, which in the time domain is
% boxcar multiplication, and shifting the boxcar into the time range
% of 0->tmax.
%
% f    ... vector of frequency coordinates. Must be regularly sampled and
%	   in increasing order.
% v    ... vector of complex spectral values
% fint ... vector of interpolation frequencies. All fint must lie between
%	   f(1) and f(length(f).
% vint ... vector of complex spectral values interpolated linearly at the
% 	   frequencies fint.
%
% G.F. Margrave and J. Bancroft, CREWES Project, U of Calgary, 1996
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

df=f(2)-f(1);
i1 = floor( fint/df) +1;
int = fint/df+1;

del = (int-i1);
v1 = v(i1).*exp(-i*pi*(del));
v2=v(i1+1).*exp(i*pi*(1-del));

vint = v1 + del.*(v2-v1);

