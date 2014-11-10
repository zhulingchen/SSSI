function [trout,x]= deconw(trin,trdsign,n,stab)
% [trout,x]=deconw(trin,trdsign,n,stab)
% [trout,x]=deconw(trin,trdsign,n)
%
% routine performs a Weiner style deconvolution of the
% input trace
%
% trin= input trace to be deconvolved
% trdsign= input trace to be used for operator design
% n= number of autocorrelogram lags to use (and length of
%    inverse operator
% stab= stabilization factor expressed as a fraction of the
%       zero lag of the autocorrelation.
%      ********* default= .0001 **********
%
% trout= output trace which is the deconvolution of trin
% x= output inverse operator used to deconvolve trin
%
% by: G.F. Margrave, May 1991
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
  if nargin<4
    stab=.0001;
  end
% generate the autocorrelation
  a=auto(trdsign,n,0);
% stabilize the auto
  a(1)=a(1)*(1.0 +stab);
  a=a/a(1);
% generate the right hand side of the normal equations
  b=[1.0 zeros(1,length(a)-1)];
% do the levinson recursion
  x=levrec(a,b);
% normalize the inverse operator
  x=x/sqrt(x'*x);
% deconvolve trin
  trout=convm(trin,x);
  trout=balans(trout,trin);
  
