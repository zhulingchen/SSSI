function prfilt=predict(trin,nop,nlag,stab)

% prfilt= predict(trin,nop,nlag,stab)
% prfilt= predict(trin,nop,nlag)
% prfilt= predict(trin,nop)
%
% PREDICT returns the nop long Wiener prediction filter for a  
% lag of nlag time samples. The filter is designed on the input
% trace, trin, and a stab factor is included by multiplying the 
% zero lag of the normalized autocorrelation by 1+stab.
%
% trin= input trace used to design the prediction operator
% nop= number of points in the prediction operator
% nlag= prediction lag distance in samples
% *********** default= 1 **************
% stab= stabilazation factor expressed as a fraction of the zero
%       lag of the autocorrelation.
%  ************ default= .0001 ***********
% prfilt= minimum phase Wiener prediction filter. Designed
%  such that:  prediction error= p.e.=> 
%      p.e.=trin(nlag:length(trin))-trinhat(1:length(trin)-nlag+1)
%  has minimum squared length and where:
%      trinhat= conv(trin,prfilt) is the predictable part of trin
%
% Re: Peacock and Treitel, Geophysics vol 34, 1968
% 
% by: G.F. Margrave, July 1991
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
  if nargin<3
     nlag=1;
  end
% generate the autocorrelation
  a=auto(trin,nlag+nop,0);
% stabilize the auto
  a(1)=a(1)*(1.0 +stab);
  a=a/a(1);
% generate the right hand side of the normal equations
  b=a(nlag+1:nlag+nop);
% do the levinson recursion
  prfilt=levrec(a(1:nop),b);
  

