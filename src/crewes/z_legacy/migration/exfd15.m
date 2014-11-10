function [arymig,tmig,xmig]=exfd151(aryin,aryvel,t,x,tau)
% [arymig,tmig,xmig]=psmig1(aryin,aryvel,t,x,tau)
%
% PSMIG is a phase shift time migration routine.
%
% aryin ... matrix of zero offset data. One trace per column.
% aryvel ... velocity information. The are 2 possibilities:
%		1) if a scalar, then a constant velocity migration with
%		velocity=aryvel is performed.
%		2) if a vector, then it must be the same length as the number
%		of rows in aryin. In this case it is assumed to be an rms 
%		velocity function (of time) which is applied at all positions
%		along the section.
% t ... if a scalar, this is the time sample rate in SECONDS.
%		If a vector, it gives the time coordinates for the rows of 
%		aryin.
% x ... if a scalar, this is the spatial sample rate (in units 
%		consistent with the velocity information. If a vector, then
%		it gives the x coordinates of the columns of aryin
% tau ... a scalar indicate the step length in time, unit is millisecond.
%
%
% OUTPUT arguments
%
% arymig ... the output migrated time section
% tmig ... t coordinates of migrated data
% xmig ... x coordinates of migrated data
%
% By Xinxiang Li,  CREWES Project, U of Calgary, 1996
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
%flops(0);
[nsamp,ntr]=size(aryin);
% check the validity input arguments
%  ---- check t  ----
if(length(t)>1)
	if(length(t)~=nsamp)
		error('Incorrect time specification')
	end
	[nrow,nvol] = size(t) ;
	if nrow < nvol
		t = t' ;
	end
	dt=t(2)-t(1);
else
	dt=t;
	t=((0:nsamp-1)*dt)';
end
%  ---- checck x ----
if(length(x)>1)
	if(length(x)~=ntr)
		error('Incorrect x specification')
	end
	[nrow,nvol] = size(x) ;
	if nrow > nvol
		x = x' ;
	end
	dx=x(2)-x(1);
else
	dx = x;
	x=(0:ntr-1)*dx;
end
tmig = t;
xmig = x;
%  ---- test velocity info ----
[nvsamp,nvtr]=size(aryvel);
if(nvsamp==1 & nvtr~=1)
	%might be transposed vector
	if(nvtr==nsamp)
		aryvel=aryvel';
	else
		error('Velocity vector is wrong size');
	end
	%make velocity matrix
	aryvel=aryvel*ones(1,ntr);
elseif( nvsamp==1 & nvtr==1)
	aryvel=aryvel*ones(nsamp,ntr);
else
	if(nvsamp~=nsamp)
		error('Velocity matrix has wrong number of rows');
	elseif(nvtr==1)
		for j = 1:ntr
			vel(:,j) =aryvel(:);
		end
		aryvel =vel;
	else (nvtr ~= ntr)
		error('Velocity matrix has wrong number of columns');
	end
end
%
%ok, we now have a velocity vector the same length of the input trace
%
%  ---- check tau ----
vmax = max(max(aryvel));
%if tau >= 8.*dx^2/vmax^2/dt
%	error( ' the step tau should be less ');
%end
%arymig = zeros(size(tmig*xmig));
arymig = aryin;
clear aryin;
ktausm = ceil(tau/dt);
ndown = nsamp/ktausm;
% **** LOOP OVER downword steps ****
for idown = 1:ndown
%disp (['The step # ' int2str(idown) ' in ' int2str(ndown)]);
	top1 = (idown - 1) * ktausm;
	top2 = idown * ktausm;
	if (top2 > nsamp-1) top2 = nsamp-1; end
	aryn = zeros(2,ntr);
	aryo = zeros(2,ntr);
	
	for tlev = nsamp-1:-1:top1
		taueff = tau;
		if (tlev < top2 )
			taueff = tau*(tlev - top1)/ktausm;
		end
		aryo(1,:) = arymig(tlev+1,:);
		vel = 0.5*aryvel((top1+top2)/2+1, :);
		w = vel.*vel*taueff* dt / (4.*dx^2);
		new = aryn(2,1:ntr-2) - 2.*aryn(2,2:ntr-1) + aryn(2,3:ntr);
		old = aryo(1,1:ntr-2) - 2.*aryo(1,2:ntr-1) + aryo(1,3:ntr);
		aryn(1,2:ntr-1) = w(2:ntr-1).*(new+old)+aryo(1,2:ntr-1)+aryn(2,2:ntr-1)-aryo(2,2:ntr-1);
		
		% two end points
		aryn(1,1) = aryn(1,2);
		aryn(1,ntr)=aryn(1,ntr-1);
		
		% go to next level
		aryn(2,:)=aryn(1,:);
		aryo(2,:)=aryo(1,:);
		
		%*******************
		arymig(tlev+1,:)=aryn(1,:);
	end
	
end
		
%disp(['Total floating operation --' int2str(flops)]);
