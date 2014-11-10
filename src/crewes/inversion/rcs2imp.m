function imp=rcs2imp(rcs,znot,alg)
% imp=rcs2imp(rcs,znot)
%
% RCS2IMP takes a reflectivity vector specified in time and computes
% a corresponding impedence by integrating and exponentiating it. The
% underlying relation is the approximation that 
%  rcs ~ .5 d(log(imp))=> imp=znot*exp(2*Integral(rcs))
%
% rcs ... input column vector of reflection coefficients
% znot ... scalar value of first impedence
% imp ... the resulting impedance
%
% G.F. Margrave
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
    alg=1;
end
 %integrate the rcs
 if(alg>1)
     ircs=zeros(size(rcs));
     %generate vector of the average of adjacent points in rcs
     n=length(rcs);
     rcs2= (rcs(2:n)+rcs(1:n-1))/2.;
     ircs(2:n)=cumsum(rcs2);
 else
     ircs=cumsum(rcs);%simpler than above
 end
 %exponentiate
 imp=znot*exp(2*ircs);
