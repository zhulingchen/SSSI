function ylogout=log_datumadj(ylog,xlog,xd,yd,xdold,ydold)
% ylogout=log_datumadj(ylog,xlog,xd,yd,xdold,ydold)
% ylogout=log_datumadj(ylog,xlog,xd,yd)
%
% LOG_DATUMADJ performs datum adjustment on a log's vertical coordinates.
% The datum shift is accomplished
% by removing the old datum and installing the new one. The
% installation of a datum is done
% by subtracting the y coordinate of the datum at the logs x coordinate
% from the log's y coordinate vector.
% Datum removal is the opposite of datum installation.
% Datum shifts are rounded to the nearest whole sample (sample rate is 
%	taken to be ylog(2)-ylog(1) )
%
% ylog ... vector of y (vertical) coordinates for log
% xlog ... scalar giving the logs x coordinate
% (xd,yd) ... piecewise linear specification of the new datum
% (xdold,ydold) ... piecewise linear specification of the old
%		datum
%
% Default for the old datum is all zeros
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
if(nargin< 5 )
	xdold=xd;
	ydold=zeros(size(xd));
end
	%compute the shifts
	%first the shifts to remove the old datum
	dy=ylog(2)-ylog(1);
	delyold=interpextrap(xdold,ydold,xlog,0);
	delyold=dy*round(delyold/dy);
	%now to install the new one
	dely=interpextrap(xd,yd,xlog,0);
	dely= dy*round( dely/dy );
	%combine the shifts
	dely=dely-delyold;
	ylogout=ylog-dely;
