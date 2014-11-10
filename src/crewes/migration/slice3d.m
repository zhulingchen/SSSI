function tslice=slice3d(seis,t,y,x,times)
% SLICE3D ... extract a time slice, using interpolation, from a 3D volume
%
% tslice=slice3d(seis,t,y,x,times)
%
% seis ... 3d volume in the form t by y by x
% t ... time coordinates of seis (length(t)=size(seis,1))
% y ... y coordinates of seis (length(y)=size(seis,2))
% x ... x coordinates of seis (length(x)=size(seis,3))
% times ... time to slice at. May be a scalar in which case the slice is at
%      a constant time. If not a scalar, then it must be a matrix of size
%      (length(y),length(x)) giving the slice time at each (y,x).
% tslice ... matrix of size (y,x) of the slice
%      Entries are computed by sinc function interpolation in time.
%

[nt,ny,nx]=size(seis);
if( (nt-1)*(nx-1)*(ny-1) == 0)
    error('seismic matrix is not 3D')
end
if( length(t)~=nt)
    error(' time vector is the wrong size')
end
if( length(y)~=ny)
    error(' y vector is the wrong size')
end
if( length(x)~=nx)
    error(' x vector is the wrong size')
end
if(length(times)~=1)
    if(size(times,1)~=ny || size(times,2)~=nx)
        error('times argument is the wrong size')
    end
end

tslice=zeros(ny,nx);

for kx=1:nx
    for ky=1:ny
        if(length(times)==1)
            targettime=times;
        else
            targettime=times(ky,kx);
        end
        tslice(ky,kx)=sinci(seis(:,ky,kx),t,targettime);
    end
end