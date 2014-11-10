function [trout,tout]=strechlog(trin,t,zonein,zoneout)
% [trout,tout]=strechlog(trin,t,zonein,zoneout)
%
% STRECHLOG stretches a single time zone on a log into another. This done by
% sinc function resampling with zero phase anti-alias filtering if needed.
% trin = input trace or log (regularly sampled)
% t= vertical coordinate vector for trin. Ostensibly time but could equally well be
%		depth
%  ************** trin and t must be the same length ************
% zonein = time zone boundaries for the time zone to be streched. This should
%		be a two element vector giving [tmin tmax]
% zoneout = the input log, over the time zone zonein, will be streched to
% 			fill this time zone at the same sample rate as input
% The first output sample will be at zoneout(1)
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
% convert to column vectors
[l,m]=size(trin);
trin=trin(:);
t=t(:);
	
tmin1=zonein(1);
tmin2=zoneout(1);
tmax1=zonein(2);
tmax2=zoneout(2);
dt=t(2)-t(1);
	% compute number of output samples
	nsampout=round((tmax2-tmin2)/dt)+1;
		
	% determine the effective output sample rate
		
	dtout= (tmax1-tmin1)/(nsampout-1);
		
	trout=resamp(trin,t,dtout,zonein,0);% use zero phase anti-alias
	tout=tmin2:dt:tmax2;
		
	if(l==1)
		trout=trout';
	tout=tout';
 end
