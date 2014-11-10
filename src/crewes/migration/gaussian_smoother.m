function velsmo=gaussian_smooother(vel,x,z,hw)
% GAUSSIAN_SMOOTHER: smooth a velocity model by conv2 with a Gaussian
%
% velsmo=gaussian_smooother(vel,x,z,hw)
%
% Function smooths a velocity model by convolving it with a 2D Gaussian.
% The smoothing is actually done by convolving the Gaussian with slowness
% (1./vel) and then inverting the result. This has the effect of roughly
% conserving traveltime.
% The Gaussian is specified with a certain halfwidth and is truncated at
% two halfwidth. The velocity model is extended for two halfwidth's in
% every direction (constant extrapolation) to avoid edge effects. After
% convolution it is truncated to the original size
% vel ... 2D velocity matrix
% x ... vector of x coordinates for the velocity model. Length must be the
%       same as the number of columns of vel. Alternativly, this may be
%       specified as a single scalar in which case it is the horizontal
%       grid size
% z ... similar to x but for depth coordinates
% hw ... half-width of the Gaussian smoother in consistent spatial units
% velsmo ... smoothed velocity matrix the same size as vel
%

[nz,nx]=size(vel);
if(length(x)==1)
    dx=x;
    x=(0:nx-1)*dx;
else
    if(length(x)~=nx)
        error('x is the wrong size')
    end
    x=x(:)';
    dx=x(2)-x(1);
end
if(length(z)==1)
    dz=z;
    z=(0:nz-1)*dz;
    z=z';
else
    if(length(z)~=nz)
        error('z is the wrong size')
    end
    z=z(:);
    dz=z(2)-z(1);
end

%make a 2D Gaussian the hw wide
ngx2=round(3*hw/dx);
xg=(-ngx2:ngx2)*dx;
ngz2=round(3*hw/dz);
zg=((-ngx2:ngx2)*dx)';
sq_dist_from_cntr=(ones(size(zg))*xg).^2 + (zg*ones(size(xg))).^2;
gauss=exp(-sq_dist_from_cntr/(hw*hw));
% extend the velocity model
velextend=[vel(:,1)*ones(1,ngx2) vel vel(:,end)*ones(1,ngx2)];
velextend=[ones(ngz2,1)*velextend(1,:);velextend;ones(ngz2,1)*velextend(end,:)];
%
% now smooth
temp=conv2(1./velextend,gauss,'same')/sum(gauss(:));
velsmo=1./temp(ngz2+1:ngz2+length(z),ngx2+1:ngx2+length(x));