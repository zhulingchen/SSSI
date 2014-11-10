function ind=newest(t1,t2)
% ind=newest(t1,t2)
% 
% NEWEST compares two times (from clock or fix(clock)) and determines
% which is the most recent. If t1, then 1 is returned while a 2 
% indicates that t2 is more recent. If the times are identical, then
% 0 is returned.
%
% G.F. Margrave, April 1994
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
% compare years
if( t1(1)~=t2(1) )
	if(t1(1)<t2(1))
		ind=2;
	else
		ind=1;
	end
	return;
end
%compare months
if( t1(2)~=t2(2))
	if(t1(2)<t2(2))
		ind=2;
	else
		ind=1;
	end
	return;
end
delt= etime(t1,t2);
if( delt>0 )
	ind=1;
elseif(delt<0)
	ind=2;
else
	ind=0;
end
