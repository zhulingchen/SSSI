function vint=vrms2vint(vrms,t,flag)
% VRMS2VINT: convert rms to interval velocity
%
% vint=vrms2vint(vrms,t,flag)
% vint=vrms2vint(vrms,t)
%
% flag=0 ... return nonphysical interval velocities as NaN
%     =1 ... interpolate interval velocities from neighbors to
%            replace non-physical results
%  ******* default = 0 **********
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

if(nargin<3)
	flag=0;
end

%force column vectors
vrms=vrms(:);
t=t(:);

%compute interval velocity squared
vint=zeros(size(vrms));
nt=length(t);
i1=1:nt-1;i2=2:nt;
vrms2=vrms.^2;
vint(i1)= (vrms2(i2).*(t(i2)-t(1))-vrms2(i1).*(t(i1)-t(1)))./(t(i2)-t(i1));
%find and process non-physical ones
ind=find(vint<0);
if(~isempty(ind))
	if(flag)
		ilive=find(vint>0);
		vint(ind)=interpextrap(t(ilive),vint(ilive),t(ind));
	else
		vint(ind)=nan*ones(size(ind));
	end
end
%compute interval velocity
vint=sqrt(vint);
vint(nt)=vint(nt-1);
