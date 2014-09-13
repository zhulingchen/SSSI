function [model, snapshot] = rvsTimeSpongeFor2dAw(v, data, nBoundary, dz, dx, dt)
%
% RVSTIMESPONGEFOR2DAW Simulate 2-d acoustic wave reverse propagation using
% finite difference in time domain with the following partial differential
% equation (PDE)
%
% (1/v^2)*(d^2)u(z, x, t)/dt^2 + f(z, x, t) = (d^2)u(z, x, t)/dz^2 + (d^2)u(z, x, t)/dx^2
%
% with absorbing boundary condition (ABC): Sponge
%
% input arguments
% v(nz,nx)          velocity model
% data(nz,nx)       received data on the surface
% nBoundary         thickness of the absorbing boundary
% dx                horizontal distance per sample
% dz                depth distance per sample
% dt                time difference per sample
%
% output arguments
% model             final pressure wavefield u(z, x, 1) in reverse time domain
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
[nz, nx] = size(v);
[~, nt] = size(data);
rtm = zeros(nz+2, nx+2, 3);  % rtm(:, :, 1) - past; rtm(:, :, 2) - present; rtm(:, :, 3) - future


%% Absorbing boundary condition (ABC): Sponge
iz = 0:nBoundary-1;
% alpha = 0.015;
alpha = 0.015;
boundary = (exp(-( (alpha*(nBoundary-iz)).^2 ) )).^10;
boundary = boundary.';


%% Reverse-Time Migration

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

% filenameRtm2d = './videos/rtm2d.avi';
% objVideoRtm2d = VideoWriter(filenameRtm2d);
% open(objVideoRtm2d);

for it = nt:-1:1
    
    % update the source map for the next time step
    source = zeros(nz, nx);
    source(1, :) = data(:, it).';
    
    % finite differencing after zero-padding
    
    rtm(izi, ixi, 3) = b .* rtm(izi, ixi, 2) - rtm(izi, ixi, 1) + ...
        ax .* (rtm(izi, ixi+1, 2) + rtm(izi, ixi-1, 2)) + ...
        az .* (rtm(izi+1, ixi, 2) + rtm(izi-1, ixi, 2)) ...
        - vdtSq .* source;    
    
    % update fdm for next time iteration
    rtm(:, :, 1) = rtm(:, :, 2);
    rtm(:, :, 2) = rtm(:, :, 3);
    
    % apply absorbing boundary conditions to 3 sides using sponge ABC
    for ixb = (1:nBoundary)+1
        rtm(izb, ixb, 1) = boundary(ixb-1) .* rtm(izb, ixb, 1);
        rtm(izb, ixb, 2) = boundary(ixb-1) .* rtm(izb, ixb, 2);
        ixb2 = nx - nBoundary + ixb;
        rtm(izb, ixb2, 1) = boundary(nx-ixb2+2) .* rtm(izb, ixb2, 1);
        rtm(izb, ixb2, 2) = boundary(nx-ixb2+2) .* rtm(izb, ixb2, 2);
        izb2 = nz - nBoundary + ixb;
        rtm(izb2, :, 1) = boundary(nz-izb2+2) .* rtm(izb2, :, 1);
        rtm(izb2, :, 2) = boundary(nz-izb2+2) .* rtm(izb2, :, 2);
    end

    % update snapshot
    snapshot(:, :, it) = rtm(izi, ixi, 2);
    
%     figure(3), imagesc(snapshot(:,:,it))
%     title(['Iteration: ',num2str(it)])
%     colorbar;
%     drawnow
%     writeVideo(objVideoRtm2d, getframe(gcf));
    
end % time loop

% close(objVideoRtm2d);

% write out final wavefield
model = rtm(:,:,1);