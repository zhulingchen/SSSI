function index=oncurve(x,y,xnot,ynot)
% index=oncurve(x,y,xnot,ynot)
% 
% ONCURVE returns zero if the point (xnot,ynot) is not determined to lie on
% the piecewise linear curve represented by the vectors x & y. If
% (xnot,ynot) does lie on the curve then the returned value is the index
% such that (xnot,ynot) lies between 
% (x(index),y(index)) & (x(index+1),y(index+1)). 
% (If (xnot,ynot) is exactly equal to a point in (x,y) then index gives
% that point.
%
% G.F. Margrave, December 1993
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
it=[];
%t=clock;
n=length(x);
index=0; % assume false
closenough=10*eps;
% test for equality
ind= find(x==xnot);
if(~isempty(ind))
	ind2=find(y(ind)==ynot);
	if(~isempty(ind2))
			index=ind(ind2);
			return;
	end
end
% find the points which bracket xnot,ynot
indx = surround(x,xnot);
indy = surround(y,ynot);
%check for a vertical segment
if isempty(x(indy+1))
	if isempty(x(indy))
		indx=indy;
	end
elseif ~isempty(x(indy+1))
	if x(indy)==x(indy+1)
		indx=indy;
	end
end
for k=1:length(indx)
	it=find(indx(k)==indy);
	if(~isempty(it))
		itest=indx(k);
		break;
	end
end
if(isempty(it))
	return;
end
	
%compute perpedicular distance
if( x(itest)~=x(itest+1) )
	m=(y(itest+1)-y(itest))./(x(itest+1)-x(itest));
	b=y(itest)-m.*x(itest);
	dtest=abs(m*xnot-ynot+b)./sqrt(m.*m+1);
else
 dtest=abs(xnot-x(itest));
end
d=min(dtest);
	
if(d<closenough)
		it=find(d==dtest);
		index=itest(it);
		%etime(clock,t)
		return;
end
	
	
