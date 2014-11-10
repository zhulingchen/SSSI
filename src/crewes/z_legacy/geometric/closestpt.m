function [xpt,ypt,inode,d]=closestpt(x,y,xnot,ynot)
% [xpt,ypt,inode,d]=closestpt(x,y,xnot,ynot)
% 
% Given a piecewise linear curve whose nodes are defined by the vectors x & y
% CLOSESTPT returns the x,y coordinates of the point on the curve which is
% closest to (xnot,ynot) (which is not necessarily on the curve). This
% algorithm works well if (xnot,ynot) is close to the curve, however, if
% it is far away, then the returned point is just the closest node.
% xpt,ypt= coordinates of the closest point
% inode = closest point either is (x(inode),y(inode)) or occurs between 
% 		(x(inode),y(inode)) & (x(inode+1),y(inode+1))
% d = euclidean distance from (xnot,ynot) to (xpt,ypt)   
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
xhat1=[]; xhat2=[]; yhat1=[]; yhat2=[];
if length(x) ~= length(y)
	error(' Vectors x and y must have same length');
end
%t=clock;
n=length(x);
% test for equality
ind= find(x==xnot);
if(~isempty(ind))
	ind2=find(y(ind)==ynot);
	if(~isempty(ind2))
			xpt=x(ind(ind2));
			ypt=y(ind(ind2));
			inode=ind(ind2);
			d=0.0;
			return;
	end
end

% find the points which bracket xnot,ynot
indx = surround(x,xnot);
indy = surround(y,ynot);
it=[];
for k=1:length(indx)
	if isempty(indy)
		it=find(isempty(indx));
	else 
		it=find(indx(k)==indy);
	end
	if(~isempty(it))
		itest=indx(k);
		break;
	end
end
if(isempty(it))
	% if no points bracket (xnot,ynot) then it is way exterior to the curve
	% we have no choice but to compute the distances to the various nodes and find 
	% the minimum
	dd=sqrt( (x-xnot).^2 + (y-ynot).^2 );
	dmin=min(dd);
    indd=find(dd==dmin);
    inode=indd;
    d=dmin;
    xpt=x(inode);
    ypt=y(inode);
    return;
    
%    
end
	
%compute perpedicular distance
if( x(itest+1)~=x(itest) )
	m=(y(itest+1)-y(itest))./(x(itest+1)-x(itest));
	b=y(itest)-m.*x(itest);
	dtest=abs(m*xnot-ynot+b)./sqrt(m.*m+1);
	d=min(dtest);
	it=find(d==dtest);
	ind=itest(it);
	xpt=(m*(ynot-y(ind))+xnot+m*m*x(ind))/(m*m+1);
	ypt=m*(xpt-x(ind))+y(ind);
	inode=ind;
else
	xpt=x(itest);
	ypt=ynot;
	d=abs(x(itest)-xnot);
	inode=itest(1);
end
