function a=xy2arclen(x,y,xnot,ynot)
% a=xy2arclen(x,y,xnot,ynot)
% a=xy2arclen(x,y,xnot)
% a=xy2arclen(x,y)
%
% Given a piecewise linear curve whose nodes are specified by the vectors x & y
% XY2ARCLEN, in the form a=xy2arclen(x,y,xnot,ynot), computes the arclength 
% (or inline distance) from (x(1),y(1)) to each of the points given by the 
% vectors (xnot,ynot). If (xnot,ynot) does not lie on the curve then a NAN 
% is returned. In the form, a=xy2arclen(x,y), the arclength to each of the
% verticies of (x,y) is returned. Lastly, in the form a=xy2arclen(x,y,xnot)
% the arclength to each point on the curve with x coordinate xnot is
% returned.  In this form, xnot must be a single scalar value and the length
% of a will be greater than 1 if the curve is multi-valued. The 2 argument
% form is the most efficient since no validity checking must be done to
% determine if points lie on the curve or not.
%
% G.F. Margrave December 1993
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
if( nargin == 2)
	n=length(x);
	xx=x(2:n);
	yy=y(2:n);
	d=sqrt( (x(1:n-1)-xx).^2 + (y(1:n-1)-yy).^2 );
	ind=find(isnan(d));
	if(length(ind)>0)
		d(ind)=zeros(size(ind));
	end
	a=zeros(size(x));
	a(2:n)=cumsum(d);
	return;
end
if( nargin == 3 )
	if(length(xnot)>1)
		error('xnot must be a scalar when three arguments are supplied');
	end
	ynot=ycurve(x,y,xnot);
	if(length(ynot)>1)
		xnot=xnot*ones(size(ynot));
	end
end
n=length(xnot);
anodes=xy2arclen(x,y);
a=zeros(size(xnot));
for k=1:n
	m=oncurve(x,y,xnot(k),ynot(k));
	
	if(m==0)
		a(k)=NaN;
	else
	 % note the explicit use of m(1) here. If m has more than one value, then
	 % the curve is multivalued and that will not be handled correctly except int
	 % the case of logsec segments which are linked at the multi valued points
		ainc=sqrt( (x(m(1))-xnot(k)).^2 + (y(m(1))-ynot(k)).^2);
		a(k)=anodes(m(1))+ainc;
	end
	
end
	
