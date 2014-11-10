function [seis,t]=afd_shotrec_alt(dx,dtstep,dt,tmax, ...
         velocity,xshot,zshot,xrec,zrec,filt,phase,laplacian)
% AFD_SHOTREC_ALT ... makes finite difference shot records
%
% [seis,t]=afd_shotrec_alt(dx,dtstep,dt,tmax,velocity,xshot,zshot,xrec,zrec,filt,phase,laplacian)
%
%
% AFD_SHOTREC will create a shot record given a velocity model and source 
% and receiver configutations. 
% Two input snapshots of the wavefield, one at time=0-dtstep and one at 
% time=0, are used in a finite difference algorithm to propogate the 
% wavefield.  The finite difference algorithm can be calculated with a 
% five or nine point approximation to the Laplacian operator.  The five 
% point approximation is faster, but the nine point results in a 
% broader bandwidth. Note that the velocity and grid spacing must 
% fulfill the equation max(velocity)*dtstep/dx < 0.7 for the model 
% to be stable.  The source array is included in the input 
% wavefield matrices, snap1 and snap2, and the receiver locations are defined
% by the user. Receiver arrays are not explicitely provided though they
% can be simulated after-the-fact by appropriate spatial convolutions.
% Of the two required snapshots, commonly they are either set equal to each other
% or the earlier one (snap1) is set to zero. These give different initial
% conditions. Two seismograms are returned:  one including all frequencies, and
% the other filtered by a gaussian function specified by the user.
%
% dx = the bin spacing for both horizontal and vertical (in consistent units)
% dtstep = size of time step for modelling (in seconds)
%         Stability requires max(velocity)*dtstep/dx < 0.7.
% dt = size of time sample rate for output seismogram. dt<dtstep causes resampling.
%		 dt>dtstep is not allowed. This allows the model to be oversampled for
%		 propagation but then conventiently resampled.
%		 The sign of dt controls the phase of the antialias resampling filter.
%		 dt>0 gives minimum phase, dt<0 is zero phase. Resampling is of course
%		 done at abs(dt). 
% tmax = the maximum time of the seismograms in seconds
% velocity = the input velocity matrix in consisnent units
%          must have the same size as snap1 and snap2
% xshot = array of x coordinates of the desired shot points (in consisent units)
% zshot = array of z coordinates of the desired shot points (in consisent units)
% xrec = a vector of the x-positions of receivers (in consisent units)
% zrec = a vector of the z-positions of receivers (in consistent units)
% filt or wlet = a 4 component ' Ormsby ' specification = [f1 f2 f3 f4] in Hz
%        or a wavelet. If it is longer than four elements, it is assumed to be
%      a wavelet. The wavelet should be sampled at the time step sample rate
%      (dtstep)
% phase or tw ... If a the previous input was a four point filter, then this must
%     a scalar where 0 indicates a zero phase filter and 1 is a minimum phase
%     filter. If the previous input was a wavelet, then this is the time 
%     coordinate vector for the wavelet. The time sample rate of the wavelet MUST
%		equal dtstep.
% laplacian - an option between two approximation to the laplacian operator
%           - 1 is a 5 point approximation
%           STABILITY CONDITION: max(velocity)*dtstep/dx MUST BE < sqrt(2)
%           - 2 is a nine point approximation
%           STABILITY CONDITION: max(velocity)*dtstep/dx MUST BE < sqrt(3/8)
%  ************** default = 1***********
%
% seis = the output seismogram
% t = the time vector (from 0 to tmax)
%
% by Carrie Youzwishen, February 1999
%    G. F. Margrave July 2000
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

tic;
boundary=2;

[nz,nx]=size(velocity);

xmax=(nx-1)*dx;
zmax=(nz-1)*dx;

x=0:dx:xmax;
z=(0:dx:zmax)';

nrec=length(xrec);
if(nrec~=length(zrec))
	error('xrec and zrec are inconsistent')
end

nshot=length(xshot);
if(nshot~=length(zshot))
	error('xshot and zshot are inconsistent')
end

test=between(0,xmax,xrec,2);
if(length(test)~=length(xrec))
	error('xrec not within grid')
end

test=between(0,zmax,zrec,2);
if(length(test)~=length(zrec))
	error('zrec not within grid')
end
	

if laplacian ==1 
    if max(max(velocity))*dtstep/dx > 1/sqrt(2)
    	error('Model is unstable:  max(velocity)*dtstep/dx MUST BE < 1/sqrt(2)');
    end
elseif laplacian ==2
    if max(max(velocity))*dtstep/dx > sqrt(3/8)
    	error('Model is unstable:  max(velocity)*dtstep/dx MUST BE < sqrt(3/8)');
    end
else
   error('invalid Laplacian flag')
end


if(abs(dt)<dtstep)
	error('abs(dt) cannot be less than dtstep')
end

if(length(filt)==4 & length(phase)==1 )
   %build Ormsby wavelet
   tlen= 8/(filt(2)+filt(3));
   [w,tw]=ormsby(filt(1), filt(2), filt(3), filt(4), tlen, dtstep);
   if(phase~=0 & phase ~=1)
      error('invalid phase flag');
  elseif(phase==1)
      w=tomin(w,.0001);
  end
else
   if(length(filt)~=length(phase))
      error('invalid wavelet specification')
   end
   w=filt;
   tw=phase;
   filt=0;
   if abs( tw(2)-tw(1)-abs(dtstep)) >= 0.000000001
      error('the temporal sampling rate of the wavelet and dtstep MUST be the same');
   end
end

nw=length(w);

%temporal information chosen by user
%t=(0:dtstep:tmax)';

%set up matrix for output seismogram
seis=zeros(floor(tmax/dtstep),nrec);

%transform receiver locations to bin locations
ixrec = floor(xrec./dx)+1;
izrec = floor(zrec./dx)+1;

%determine linear addresses for receivers
irec=(ixrec-1)*nz + izrec;

%transform source locations to bin locations
ixshot = floor(xshot./dx)+1;
izshot = floor(zshot./dx)+1;

%determine linear addresses for shots
ishot=(ixshot-1)*nz + izshot;

%determine time shift to allow for non-causal wavelets
inot=find(tw==0);
if(isempty(inot))
    error('wavelet does not hav a sample at time 0');
end

%build snapshots
snap1=zeros(size(velocity));
snap2=snap1;
snap2(ishot)=w(1);

%grab time zero from snap2
seis(1,:)=snap2(irec);

maxstep=round(tmax/dtstep)-1+(inot-1);
disp(['There are ' int2str(maxstep) ' steps to complete']);
time0=clock;

% each loop does two time steps
nwrite=2*round(maxstep/50)+1;
for k=1:2:maxstep
	
	%time step
	snap1=afd_snap(dx,dtstep,velocity,snap1,snap2,laplacian,boundary);
	seis(k+1,:)=snap1(irec);
    if((k+1)<=nw)
        snap2(ishot)=w(k+1);
    end
    
    if(iscomplex(snap1))
        yyy=1;
    end
	
	snap2=afd_snap(dx,dtstep,velocity,snap2,snap1,laplacian,boundary);
	seis(k+2,:)=snap2(irec);
    if((k+2)<=nw)
        snap2(ishot)=w(k+2);
    end
    

        if rem(k,nwrite) == 0
           timenow=clock;
           tottime=etime(timenow,time0);

           disp(['wavefield propagated to ' num2str(k*dtstep) ...
           ' s; computation time left ' ...
            num2str((tottime*maxstep/k)-tottime) ' s']);
        end 

end

%get rid of negative time samples
if(inot>1)
    seis(1:inot-1,:)=[];
end

%compute a time axis
t=((0:size(seis,1)-1)*dtstep)';

disp('modelling completed')

%resample if desired
if(abs(dt)~=dtstep)
	disp('resampling')
	phs=(sign(dt)+1)/2;
	dt=abs(dt);
	for k=1:nrec
		cs=polyfit(t,seis(:,k),4);
		[tmp,t2]=resamp(seis(:,k)-polyval(cs,t),t,dt,[min(t) max(t)],phs);
		seis(1:length(tmp),k)=tmp+polyval(cs,t2);
	end
	seis(length(t2)+1:length(t),:)=[];
	t=t2;
end
if(iscomplex(seis))
        %disp('Really really bad! Complex amplitudes generated...')
        seis=real(seis);
end

toc;
