function arymig=kirkcv(aryin,v,dt,dx,aper)
%
% arymig=kirk(aryin,v,dt,dx,aper)
%
% KIRKCV is a constant velocity post stack Khirchoff time 
%	migration routine.
%
% aryin ... matrix of zero offset data. One trace per column.
% v ... velocity
% dt ... the time sample rate in SECONDS
% dx ... the spatial sample rate in units consistent with the velocity
%		information.
% aper ... migration aperture in physical length units
% arymig ... the output migrated time section
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
t1=clock; % save start time
[nsamp,ntr]=size(aryin);
%compute half-aperture in traces
traper=round(.5*aper/dx);
%initialize output array
arymig=zeros(nsamp,ntr);
%one way time
dt2=.5*dt;
t=((0:nsamp-1)*dt2)';
%compute maximum time needed
tmaxin=2*t(nsamp);
tmax=sqrt( tmaxin^2 + (aper/v)^2);
%pad input to tmaxin
npad=1.05*(tmax-tmaxin)/dt;
aryin= [aryin; zeros(npad,ntr)];
%precompute hyperbolic lookup table
htab=zeros(nsamp,traper+1);
htab(:,1)=(1:nsamp)';
zo2=(t*v).^2;
for k=1:traper
	offset2=(k*dx).^2;
	htab(:,k+1)= sqrt(zo2 + offset2)./(dt2*v);
end
%loop over migrated traces
for ktr=1:ntr
	%determine traces in aperture
	n1=max([1 ktr-traper]);
	n2=min([ntr ktr+traper]);
	truse=n1:n2;
	
	%offsets and depths
	offset2=((truse-ktr)*dx).^2;
	zo2=(t*v).^2;
	
	%loop over migrated samples
	for ksamp=1:nsamp
		%compute contributing sample numbers across aperture
		%ksamp_use= sqrt(zo2(ksamp) + offset2)./(dt2*v);
		
		%check for samples out of bounds
		%ind=find(ksamp_use>nsamp);
		%if(~isempty(ind))
		%	ksamp_use(ind)=zeros(size(ind));
		%end
		ksamp_use=htab(ksamp,abs(truse-ktr)+1);
		
		%sum across hyperbola
		arymig(ksamp,ktr) = sum( matraj(aryin(:,truse),ksamp_use) );
	end
end
t1=etime(clock,t1);
disp(['Total elapsed time ' num2str(t1)])
