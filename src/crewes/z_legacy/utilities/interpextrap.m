function y0=interpextrap(x,y,x0,flag)
% y0=interpextrap(x,y,x0,flag)
% y0=interpextrap(x,y,x0)
%
% INTERPEXTRAP is idential to MATLAB's INTERP1 for x0 which lie within
% the bounds of x. For x0 < min(x) or x0> max(x) a linear extrapolation
% is done using the constant slope of the closest segment of (x,y) (flag==1)
% or a horizontal extrapolation (flag==0). Default for flag is 1
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
if(nargin<4) flag=1; end
y0=nan*ones(size(x0));
nx=length(x);
%test for trivial case
if( length(x)==1 )
	y0=y;
	return;
end
% handle the end extrapolation 
if( x(1)<x(nx) )% case for a normally ordered x
	ind=find( x0< x(1) );% points in the beginning
	if(~isempty(ind))
		if(flag)
			m1=(y(2)-y(1))/(x(2)-x(1));
		else
			m1=0.;
		end
		y0(ind)=m1*(x0(ind)-x(1))+y(1);
	end
	ind=find( x0>x(nx) );%points at the end
	if(~isempty(ind))
		if(flag)
			m2=(y(nx)-y(nx-1))/(x(nx)-x(nx-1));
		else
			m2=0.;
		end
		y0(ind)=m2*(x0(ind)-x(nx))+y(nx);
	end
else
	ind=find( x0> x(1) );
	if(~isempty(ind))
		if(flag)
			m1=(y(2)-y(1))/(x(2)-x(1));
		else
			m1=0.;
		end
		y0(ind)=m1*(x0(ind)-x(1))+y(1);
	end
	ind=find( x0<x(nx) );
	if(~isempty(ind))
		if(flag)
			m2=(y(nx)-y(nx-1))/(x(nx)-x(nx-1));
		else
			m2=0.;
		end
		y0(ind)=m2*(x0(ind)-x(nx))+y(nx);
	end
end
ind=isnan(y0);
y0(ind)=lint(x,y,x0(ind));
	
