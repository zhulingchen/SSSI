function w = mwindow(n,percent)
% MWINDOW: creates an mwindow (boxcar with raised-cosine tapers)
%
% w = mwindow(n,percent)
% w = mwindow(n)
% 
% MWINDOW returns the N-point Margrave window in a 
% column vector. This window is a boxcar over the central samples
% (100-2*percent)*n/100 in number, while it has a raised cosine
% (hanning style) taper on each end. If n is a vector, it is
% the same as mwindow(length(n)
%
% n= input length of the mwindow. If a vector, length(n) is
%    used
% percent= percent taper on the ends of the window
%   ************* default=10 ************
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
  
% set defaults
 if nargin<2
  percent=10;
 end
 if length(n)>1
   n=length(n);
 end
% compute the Hanning function 
 if(percent>50)|(percent<0)
   error(' invalid percent for mwindow')
 end
 m=2.*percent*n/100.;
 m=2*floor(m/2);
 h=hanning(m);
 w = [h(1:m/2);ones(n-m,1);h(m/2:-1:1)];



