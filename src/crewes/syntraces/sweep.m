function [s,t]=sweep(fmin,fmax,dt,tmax,taper)
% SWEEP: generate a linear Vibroseis sweep
%
% [s,t]= sweep(fmin,fmax,dt,tmax,taper)
% [s,t]=sweep(fmin,fmax,dt,tmax)
%
% SWEEP generates a linear synthetic Vibroseis sweep for the 
% specified passband. (After Waters, Reflection Seismology, 1981, page 90)
%
% fmin= minimum swept frequency in Hz
% fmax= maximum swept frequency in Hz
% dt= time sample rate in seconds
% tmax= sweep length in seconds
% taper= length of cosine sweep taper in seconds
% ********** default =.5 if tmax>4.0
%            otherwise =.25 *************
%
% s= output linear sweep
% t= output time coordinate vector 
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

 if nargin<5
   if tmax>4.0
    taper=.5;
   else
    taper=.25;
   end
 end
%
 t=0.:dt:tmax;
 b=(fmax-fmin)/tmax;
% 
 f=fmin+.5*b*t;
 theta=2.*pi*f.*t;
 s=sin(theta);
% apply taper
 percent=100.*taper/tmax;
 s=s.*mwindow(s,percent)';

 
