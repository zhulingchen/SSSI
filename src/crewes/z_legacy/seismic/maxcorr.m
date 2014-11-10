function x=maxcorr(trace1,trace2,n)
% x=maxcorr(trace1,trace2,n)
%
% MAXCORR computes 2*n+1 lags of the crosscorrelation of trace1
% with trace2 (using CCORR) and then uses splines to pick the
% maximum to the nearest .1 lag. The interpolated maximum and 
% lag are returned as x(1) and x(2). Here "maximum" means maximum
% absolute value.
%
% trace1= input trace number 1
% trace2= input trace number 2
% n= 2*n +1 lags will be computed
% x= output: x(1)-> interpolated maximum cross correlation
%            x(2)-> interpolated lag at maximum correlation
% 
% Note: a negative result for x(2) indicates trace2 is delayed
%       relative to trace 1
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
a=ccorr(trace1,trace2,n);
 amax=max(abs(a));
 imax=near(abs(a),amax);
 amax=a(imax);
% interpolate a maximum with splines
 x1=max(1,imax-2);
 x2=min(length(a),imax+2);
 xs=x1:x2;
 xi=x1:.1:x2;
 ai=spline(xs,a(xs),xi);
 amax=max(abs(ai));
 imax=near(abs(ai),amax);
 amax=ai(imax);
 x=[amax xi(imax)-n-1];
