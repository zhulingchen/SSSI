function [data, snapshot] = fwdTimeSpongeFor2dAw(v, source, nBoundary, dz, dx, dt)
%
% FWDTIMESPONGEFOR2DAW Simulate 2-d acoustic wave forward propagation using
% finite difference in time domain with the following partial differential
% equation (PDE)
%
% (1/v^2)*(d^2)u(z, x, t)/dt^2 + f(z, x, t) = (d^2)u(z, x, t)/dz^2 + (d^2)u(z, x, t)/dx^2
% 
% with absorbing boundary condition (ABC): Sponge
%
% input arguments
% v(nz, nx)         velocity model
% source(nz, nx)    source vector (e.g., shots)
% nBoundary         thickness of the absorbing boundary
% dx                horizontal distance per sample
% dz                depth distance per sample
% dt                time difference per sample
%
% output arguments
% data              received 2-d x - time signal on the surface
% snapshot          pressure field u(z, x, t) in time domain
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Originally written by
% Stuart Kozola for Large Data in MATLAB: A Seismic Data Processing Case Study
% http://www.mathworks.com/matlabcentral/fileexchange/30585-large-data-in-matlab--a-seismic-data-processing-case-study
% The MathWorks Inc., 03/2011
% Modified by
% Lingchen Zhu (zhulingchen@gmail.com), Entao Liu (liuentao@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

%% Initialize storage
[nz, nx, nt] = size(source);
data = zeros(nx, nt);
fdm = zeros(nz+2, nx+2, 3);  % fdm(:, :, 1) - past; fdm(:, :, 2) - present; fdm(:, :, 3) - future


%% Absorbing boundary condition (ABC): Sponge
iz = 0:nBoundary-1;
% alpha = 0.015;
alpha = 0.015;
boundary = (exp(-( (alpha*(nBoundary-iz)).^2 ) )).^10;
boundary = boundary.';


%% Forward-Time Modeling

% finite difference coefficients
ax = (v*dt/dx).^2;
az = (v*dt/dz).^2;
vdtSq = (v*dt).^2;
b = 2-2*ax-2*az;

% common indicies
izi = 2:(nz+1);     % used for z interior of fdm after zero-padding
ixi = 2:(nx+1);     % used for x interior of fdm after zero-padding
izb  = (1:nz-nBoundary)+1;      % boundary z

% snapshot
snapshot = zeros(nz, nx, nt);

% filenameFm2d = './videos/fm2d.avi';
% objVideoFm2d = VideoWriter(filenameFm2d);
% open(objVideoFm2d);

for it = 1:nt
    
    % finite differencing after zero-padding
    
    fdm(izi, ixi, 3) = b .* fdm(izi, ixi, 2) - fdm(izi, ixi, 1) + ...
        ax .* (fdm(izi, ixi+1, 2) + fdm(izi, ixi-1, 2)) + ...
        az .* (fdm(izi+1, ixi, 2) + fdm(izi-1, ixi, 2)) ...
        - vdtSq .* source(:, :, it);    
    
    % update fdm for next time iteration
    fdm(:, :, 1) = fdm(:, :, 2);
    fdm(:, :, 2) = fdm(:, :, 3);
    
    % apply absorbing boundary conditions to 3 sides using sponge ABC
    for ixb = (1:nBoundary)+1
        fdm(izb, ixb, 1) = boundary(ixb-1) .* fdm(izb, ixb, 1);
        fdm(izb, ixb, 2) = boundary(ixb-1) .* fdm(izb, ixb, 2);
        ixb2 = nx - nBoundary + ixb;
        fdm(izb, ixb2, 1) = boundary(nx-ixb2+2) .* fdm(izb, ixb2, 1);
        fdm(izb, ixb2, 2) = boundary(nx-ixb2+2) .* fdm(izb, ixb2, 2);
        izb2 = nz - nBoundary + ixb;
        fdm(izb2, :, 1) = boundary(nz-izb2+2) .* fdm(izb2, :, 1);
        fdm(izb2, :, 2) = boundary(nz-izb2+2) .* fdm(izb2, :, 2);
    end
    
    % update data
    data(:, it) = fdm(2, ixi, 2);
    
    % update snapshot
    snapshot(:, :, it) = fdm(izi, ixi, 2);

%     figure(3), imagesc(snapshot(:,:,it))
%     title(['Iteration: ',num2str(it)])
%     colorbar;
%     drawnow
%     writeVideo(objVideoFm2d, getframe(gcf));

end % time loop

% close(objVideoFm2d);

%data = data(21:nx-nBoundary,:);