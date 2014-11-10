function s=slope(x,y,xnot)
% s=slope(x,y,xnot)
%
% compute the slope of the curve given by x,y at the point xnot. If the curve
% is multivalued, then multiple s's will be returned. If the curve is a polygon, then
% it must have its first and last points explicitly equal. The curve is assumed 
% piecewise linear. If xnot is exactly equal to a point in x, then the slope 
% returned is the average of the segments before and after xnot.
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
	ind=find(x==xnot);
	if(isempty(ind))
		ind=surround(x,xnot);
		s=(y(ind)-y(ind+1))./(x(ind)-x(ind+1));
	else
			n=length(x);
			% extrapolate ends for end conditions
			m1=(y(2)-y(1))/(x(2)-x(1));
			x0=x(1)-(x(2)-x(1));
			y0=y(1)+(x0-x(1))*m1;
			mn=(y(n)-y(n-1))/(x(n)-x(n-1));
			xn1=x(n)+(x(n)-x(n-1));
			yn1=y(n)+(xn1-x(n))*mn;
			x=[x0 x xn1];% pad first and last samples for end conditions
			y=[y0 y yn1];
			ind=ind+1;
			s1=(y(ind)-y(ind+1))./(x(ind)-x(ind+1));
			s2=(y(ind-1)-y(ind))./(x(ind-1)-x(ind));
		s=.5*(s1+s2);
	end
