function ynot=ycurve(x,y,xnot)
% ynot=ycurve(x,y,xnot)
%
% compute the y coordinates of the curve given by x,y at the point xnot.
% If the curve is multivalued, then multiple y's will be returned. If the
% curve is a polygon, then it must have its first and last points explicitly
% equal. If xnot lies outside the bounds and the curve is not a polygon,
% then the curve is extrapolated by repeating the% first or last sample as
% needed. If it is outside the bounds and the curve is a polygon, then []
% is returned.
%
% by G.F. Margrave December 1993
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
		if(isempty(ind)) % test for out of range
			if( x(1)==x(length(x)) & y(1)==y(length(y)) )
				ynot=[]; % don't extrapolate a polygon
				return;
			else
				% extrapolate at constant y either the first or last point
				xends=[x(1) x(length(x))];
				yends=[y(1) y(length(y))];
				if( xnot < xends(1) & xnot< xends(2) )
					ind=find(xends==min(xends));
					ynot=yends(ind);
					return;
				else
					ind=find(xends==max(xends));
					ynot=yends(ind);
					return;
				end
			end
		end
				
		ynot=y(ind)-(x(ind)-xnot).*(y(ind)-y(ind+1))./(x(ind)-x(ind+1));
	else
		ynot=y(ind);
	end
