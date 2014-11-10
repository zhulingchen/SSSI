function [r,t]=reflec(tmax,dt,amp,m,n)
% REFLEC: synthetic pseudo random reflectivity
%
% [r,t]=reflec(tmax,dt,amp,m,n)
% [r,t]=reflec(tmax,dt,amp,m)
% [r,t]=reflec(tmax,dt,amp)
% [r,t]=reflec(tmax,dt)
%
% REFLEC creates a psuedo random reflectivity by first 
% generating a Gaussian random noise sequence and then raising
% each element in the sequence to a small integral power. This 
% has the effect of making the sequence more spiky.
%
% tmax = record length
% dt= sample rate
% amp= maximum rc;
% ************* default=.2 *******************
% m= exponentiation power to which gaussian distribution 
%    is raised; Should be an odd power to preserve the sign 
%    of the sample.
% ***************** default= 3 ***********************
% n= random number seed; 
% ************* defaults to a random number based on the
%                system clock  ***********
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

if(nargin<3)
  amp=.2;
 end
if(nargin<4)
 m=3;
end
if(nargin<5)
  c=clock;
  n=c(6);
end
m=round(m);
t=(0.:dt:tmax)';
% matlab random number generator screws up with a negative seed
randn('seed',abs(n))
if(floor(m/2)*2==m)
    tmp=randn(size(t));
    r=(tmp.^m).*sign(tmp);
else
    r=randn(size(t)).^m;
end
r=amp*r/max(abs(r));
