function [val] = fpequal(v1,v2,factor)
%FPEQUAL compares two floating point numbers and returns a logical value if
%they are equal.  In the case where vectors are being compared a vector the
%same size as the input vectors will be retuned  
%
% val=fpequal(v1,v2,factor);
%
% Floating point numbers can often be different slightly eventhough they
% appear the same.  Therfore using this function will give a better
% approximation to if they are equal or not.
%
% Inputs:
%       v1 - first number or vector to be compared
%       v2 - second number or vector to be compared
%       factor - the sensitivity, larger numbers equal less sensitivity
%               ************Default 100***********
%
% If Using vectors they must be the same size as each other
%
%
%G.F. Margrave & H.J.E Lloyd 2009
%
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

if (nargin<3), factor=100; end;

small=factor*eps;
val=abs(v1-v2)<small;

end

