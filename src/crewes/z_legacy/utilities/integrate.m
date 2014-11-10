function a=integrate(u,x)
% a=integrate(u,x)
% a=integrate(u)
%
% INTEGRATE numerically integrates the vector u with respect to x.
% That means that x must be a vector of the same length as u giving
% the xcoordinates of u. If x is not supplied, one is generated as 
% x=1:length(u) which implies unit spacing between the elements of
% u. For wavelengths large compared to the grid spacing, INTEGRATE
% does an approximate job of undoing MATLAB's GRADIENT function to
% within an additive constant. The method is essentially Simpson's
% trapezoidal rule. See Numerical Recipies (2nd Ed.) eqn 4.1.11
%
%	example: 
%		x=0:2*pi/1000:2*pi;
%		y=cos(x);
%		yp=gradient(y,x);
%		yapprox=integrate(yp,x)+y(1);
%		plot(x,y,x,yp,x,yapprox,'r.');
%		rmserr=sqrt(sum((y-yapprox).^2)/(length(x)))
%
% G.F. Margrave 1994
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
	nu=length(u);
	[m,n]=size(u);
	if(nargin<2)
		x=1:nu;
		if(n==1)
			x=x';
		end
	end
	dx=x(2:nu)-x(1:nu-1);
	a=zeros(size(u));
	b=u(1:nu-1) + u(2:nu);
	a(2:nu)=b.*dx/2;
	a=cumsum(a);
