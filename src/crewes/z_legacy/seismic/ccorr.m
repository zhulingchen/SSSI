function a = ccorr(v,u,n,flag)
% a=ccorr(v,u,n,flag)
% a=ccorr(v,u,n)
%
% CCORR computes 2*n+1 lags of the cross correlation (normalized) 
% of the u with v. The zeroth lag is a(n+1). That is a=u x v where x is
% cross correlation.
% 
% v= input vector
% u= input vector
% n= number of lags desired 
% flag= 1.0 ... normalize )
%       anything else ... don't normalize
%       ******* default =1.0 ******
%
%
% NOTE: A two sided autocorrelation or a cross correlation can also
%       be computed with XCORR in the signal toolbox.
%
%   by G.F. Margrave, June 1991
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
 if (nargin<4)
  flag=1.0;
 end
% 
% master loop
% 
 [l,m]=size(v);
 lu=size(u,1);
 done=0;
% for row vectors
 if l==1 
   if lu==1, u=u'; else u=conj(u); end
   nzero=n+1;
   a(nzero)=v*u;
   uf=u;
   ub=u;
   for k=1:n
     uf=[0.0; uf(1:length(uf)-1)];
     ub=[ub(2:length(ub));0.0];
% forward lag
     a(nzero+k)=v*uf;
% backward lag
     a(nzero-k)=v*ub;
   end
   done=1;
 end
% for column vectors
 if m==1 
   v=v.';
   if lu==1, u=u'; else u=conj(u); end
   nzero=n+1;
   a(nzero)=v*u;
   uf=u;
   ub=u;
   for k=1:n
     uf=[0.0; uf(1:length(uf)-1)];
     ub=[ub(2:length(ub));0.0];
% forward lag
     a(nzero+k)=v*uf;
% backward lag
     a(nzero-k)=v*ub;
   end
   done=1; 
 end
 if done==0
   error(' input not a vector')
 end
% normalize
 if flag==1.0
   a=2.*a/(sum(v.^2)+sum(u.^2));
 end
			     
