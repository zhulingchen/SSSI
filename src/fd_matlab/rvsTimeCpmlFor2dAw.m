function [model, snapshot] = rvsTimeCpmlFor2dAw(v, data, nDiffOrder, nBoundary, dz, dx, dt)
%
% RVSTIMECPMLFOR2DAW Simulate 2-d acoustic wave reverse propagation using
% finite difference in time domain with the following partial differential
% equation (PDE)
%
% (1/v^2)*(d^2)u(z, x, t)/dt^2 + f(z, x, t) = zP + xP
% Update xPhi, zPhi, xA, zA, xPsi, zPsi, xP, zP and solve u(z, x, t) with Nonsplit Convolutional-PML (CPML)
%
% input arguments
% v(nz,nx)          velocity model
% data(nx,nt)       received data on the surface
% nDiffOrder        number of approximation order for differentiator operator
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
% Written by Lingchen Zhu (zhulingchen@gmail.com), Entao Liu (liuentao@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

%% Initialize storage
[nz, nx] = size(v);
[~, nt] = size(data);
coeff = dCoef(nDiffOrder, 's');
k = 2 * nDiffOrder - 1;
rtm = zeros(nz+2*k, nx+2*k, 3);  % rtm(:, :, 1) - future; rtm(:, :, 2) - present; rtm(:, :, 3) - past


%% Absorbing boundary condition (ABC): Nonsplit Convolutional-PML (CPML)
ixb = 1:nBoundary;          % index of x outside left boundary
ixb2 = nx-nBoundary+ixb;    % index of x outside right boundary
izb  = 1:nz-nBoundary;      % index of z inside down boundary
izb2 = nz-nBoundary+ixb;    % index of z outside down boundary

xDampLeft = dampPml(repmat(fliplr(ixb) * dx, nz, 1), v(:, ixb), nBoundary * dx);
xDampRight = dampPml(repmat(ixb * dx, nz, 1), v(:, ixb2), nBoundary * dx);
xDamp = [xDampLeft, zeros(nz, nx-2*nBoundary), xDampRight];
xb = exp(-xDamp * dt);

zDampDown = dampPml(repmat(ixb.' * dz, 1, nx), v(izb2, :), nBoundary * dz);
zDamp = [zeros(nz-nBoundary, nx); zDampDown];
zb = exp(-zDamp * dt);

% additional arrays for storage
zPhi = zeros(nz+2*k, nx);
xPhi = zeros(nz, nx+2*k);
zA = zeros(nz+2*k, nx);
xA = zeros(nz, nx+2*k);
zPsi = zeros(nz+k, nx);
xPsi = zeros(nz, nx+k);
zP = zeros(nz+k, nx);
xP = zeros(nz, nx+k);


%% Reverse-Time Migration

% finite difference coefficients
vdtSq = (v*dt).^2;

% common indicies
izi = (1+k):(nz+k);                         % used for z interior of rtm after zero-padding
ixi = (1+k):(nx+k);                         % used for x interior of rtm after zero-padding
izl = nDiffOrder:(nz+2*k-nDiffOrder);       % updated indices for zA
ixl = nDiffOrder:(nx+2*k-nDiffOrder);       % updated indices for xA

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
    
    zPhi(izi, :) = zb .* zPhi(izi, :) + (zb - 1) .* diffOperator(rtm(izl+1, ixi, 2), coeff, dz, 1);
    xPhi(:, ixi) = xb .* xPhi(:, ixi) + (xb - 1) .* diffOperator(rtm(izi, ixl+1, 2), coeff, dx, 2);
    zA(izl, :) = diffOperator(rtm(:, ixi, 2), coeff, dz, 1) + zPhi(izl, :);
    xA(:, ixl) = diffOperator(rtm(izi, :, 2), coeff, dx, 2) + xPhi(:, ixl);
    
    zPsi(izi, :) = zb .* zPsi(izi, :) + (zb - 1) .* diffOperator(zA(izl, :), coeff, dz, 1);
    xPsi(:, ixi) = xb .* xPsi(:, ixi) + (xb - 1) .* diffOperator(xA(:, ixl), coeff, dx, 2);
    zP(izi, :) = diffOperator(zA(izl, :), coeff, dz, 1) + zPsi(izi, :);
    xP(:, ixi) = diffOperator(xA(:, ixl), coeff, dx, 2) + xPsi(:, ixi);
    
    rtm(izi, ixi, 3) = vdtSq .* (zP(izi, :) + xP(:, ixi) - source) + 2 * rtm(izi, ixi, 2) - rtm(izi, ixi, 1);
    
    % *******************************************************************
    % FCT elimination of numerical dispersion. Needed for high frequency
    % source function and makes the algorithm very slow!!!)
    % *******************************************************************
    % rtm(:,:,3) = fctForAw(rtm(:,:,1),rtm(:,:,2),rtm(:,:,3));
    
    % update rtm for next time iteration
    rtm(:, :, 1) = rtm(:, :, 2);
    rtm(:, :, 2) = rtm(:, :, 3);
    
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
model = rtm(:, :, 1);