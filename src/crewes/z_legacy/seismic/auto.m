function a = auto(v,n,flag)
% a=auto(v,n,flag)
% a=auto(v,n)
% a=auto(v)
%
% auto computes n lags of the one sided autocorrelation of 
% the vector 'v'. The first lag, a(1), is termed the 'zeroth lag'
% 
% v= input vector
% n= number of lags desired (can be no larger than length(v)).
%    ********* default =length(v) *********
% flag= 1.0 ... normalize the 'zero lag' (first returned value)
%               to 1.0.
%        anything else ... don't normalize
%       ******* default =1.0 ******
%
% NOTE: A two sided autocorrelation or a cross correlation can be
%       computed with XCORR in the signal toolbox.
%
%   by G.F. Margrave, May 1991
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
 if (nargin==2)
  flag=1.0;
 end
 if (nargin==1)
   n=length(v);
   flag=1.0;
 end 
% 
% master loop
% 
 [l,m]=size(v);
 done=0;
% for row vectors
 if l==1 
   u=v';
   for k=1:n
     a(k)=v*u;
     v=[0.0 v(1:length(v)-1)];
   end
   done=1;
 end
% for column vectors
 if m==1 
   u=v';
   for k=1:n
     a(k)=u*v;
     u=[0.0 u(1:length(u)-1)];
   end
   done=1;
 end
 if done==0
   error(' input not a vector')
 end
% normalize
 if flag==1.0
   a=a/max(a);
 end
			     
