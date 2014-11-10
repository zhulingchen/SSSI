function y2=paratran(x,y,x1,y1,x2)
% y2=paratran(x,y,x1,y1,x2)
%
% PARATRAN does parallel transport of a point (x1,y1) to the x coordinate
% x2 by moving parallel to the piecewise linear curve {x,y}
% If x2 lies outside the bounds of (x,y) then the curve is extended as
% needed by simple constant extrapolation of the endpoints.
%
% G.F. Margrave January 1994
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
 y1p=ycurve(x,y,x1);
 y2p=ycurve(x,y,x2);
	if(length(y1p)>1 || length(y2p)>1 )
		disp('parallel transport fails if curve is multivalued');
        disp('curve has x coordinates:')
        x
        disp('curve has y coordinates:')
        y
        error('Curve is multivalued in y, you need to fix it before continuing')
end
y2=y2p+y1-y1p;
