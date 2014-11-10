function [theta,z] = threc(t,vo,c,z)
% THREC: Compute record-length-limit to scattering angle for v(z)
%
% theta = threc(t,vo,c,z)
% [theta,z]= threc(t,vo,c,nz)
%
% Compute the zero offset scattering angle limit imposed
% by recording time.
%
% t ... maximum recording time
% vo ... initial velocity
% c ... accelerator ( v(z) = vo +c*z )
% z ... vector of depths for which theta is computed
% theta ... vector of limiting scattering angles. One for each z.
%
% G.F. Margrave, CREWES, 1997
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

%test for automatic z calculation
if (length(z) ==1)
	nz=z;
	if( c~=0)
		zmax= vo*(exp(c*t/2)-1)/c;
	else
		zmax= vo*t/2;
	end
	z=linspace(0,zmax-100000*eps,nz);
end

if( c ~= 0)
	gamma = vo./(vo+c*z);
	arg=c*t/2;
	
	csth = (1 ./gamma -cosh(arg) )/sinh(arg);
	theta=180*acos(csth)/pi;
	
else
	zmax=vo*t/2;
	theta=180*acos(z/zmax)/pi;
	
end
