function rho=gardner(v,a,m,rho)
% rho=gardner(v,a,m,rho)
% rho=gardner(v,a,m)
%
% GARDNER estimates a density vector given an instantaneous velocity vector
% and values for the empirical parameters a and m via the Gardner relation
%	rho=a*v.^m (a times v to the mth power)
% If a partial density vector is available, then it may be provided as a 4th
% argument which must be expanded with NaN's to be the same length as v. (Take 
% care that the existing densities are paired with the proper v's.) Densities
% will only be computed with Garnders relation where they are missing and any
% non-NaN input densities are simply passed to output.
% 	v= vector of instantaneous p-wave velocities
%	a= scalar multplier
%		*********** default is .31 **********
%	m= scalar exponent
%		*********** default is .25 **********
%	rho= on input, a vector of known densities which must be expanded with 
%		NaN's to be the same size as v
%	rho= on output, aa vector of densities the same size as v and computed
%		via Gardner's relation where densities were missing on input.
%
% G.F. Margrave, March 1994
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
if(nargin<4)
	rho=nan*ones(size(v));
end
if(nargin<3)
	m=.25;
end
if(nargin<2)
	a=.31;
end
ind=find(isnan(rho));
rho(ind)= a*v(ind).^m;
