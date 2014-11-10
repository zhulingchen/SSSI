function arymig=kirk(aryin,aryvel,dt,dx,aper)
%
% arymig=kirk(aryin,aryvel,dt,dx,aper)
%
% KIRK is a simple post stack Khirchoff time migration routine.
%
% aryin ... matrix of zero offset data. One trace per column.
% aryvel ... velocity information. The are 3 possibilities:
%		1) if a scalar, then a constant velocity migration with
%		velocity=aryvel is performed.
%		2) if a vector, then it must be the same length as the number
%		of rows in aryin. In this case it is assumed to be an rms 
%		velocity function (of time) which is applied at all positions
%		along the section.
%		3) if a matrix, then it must be the same size as aryin. Here
%		it is assumed to give the rms velocity for each sample location.
% dt ... the time sample rate in SECONDS
% dx ... the spatial sample rate in units consistent with the velocity
%		information.
% aper ... migration aperture in physical length units
% arymig ... the output migrated time section
t1=clock; % save start time
[nsamp,ntr]=size(aryin);
[nvsamp,nvtr]=size(aryvel);
%test velocity info
vmin=min(aryvel);
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
	elseif(ntr~=nvtr)
		error('Velocity matrix has wron number of columns');
	end
end
%ok, we now have a velocity matrix the same size as the data matrix
%compute half-aperture in traces
traper=round(.5*aper/dx);
%initialize output array
arymig=zeros(nsamp,ntr);
%one way time
dt2=.5*dt;
t=((0:nsamp-1)*dt2)';%$$$
%compute maximum time needed
tmaxin=2*t(nsamp);
tmax=sqrt( tmaxin^2 + (aper/vmin)^2);
%pad input to tmaxin
npad=1.05*(tmax-tmaxin)/dt
aryin= [aryin; zeros(npad,ntr)];
%loop over migrated traces
imax=0;
for ktr=1:ntr
	%determine traces in aperture
	n1=max([1 ktr-traper]);
	n2=min([ntr ktr+traper]);
	truse=n1:n2;
	
	%offsets and depths
	offset2=((truse-ktr)*dx).^2;
	zo2=(t.*aryvel(:,ktr)).^2;
	
	%loop over migrated samples
	for ksamp=1:nsamp
		%compute contributing sample numbers across aperture
		ksamp_use= sqrt(zo2(ksamp) + offset2)./(dt2*aryvel(ksamp,ktr));
		if(max(ksamp_use)>imax) imax=ksamp_use; end
		
		%check for samples out of bounds
		%ind=find(ksamp_use>nsamp);
		%if(~isempty(ind))
		%	ksamp_use(ind)=zeros(size(ind));
		%end
		
		%sum across hyperbola
		arymig(ksamp,ktr) = sum( matraj(aryin(:,truse),ksamp_use) );
	end
	disp(['Trace ' int2str(ktr) ' migrated'])
end
imax
t1=etime(clock,t1);
disp(['Total elapsed time ' num2str(t1)])
