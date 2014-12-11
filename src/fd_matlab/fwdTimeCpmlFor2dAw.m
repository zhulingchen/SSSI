function [data, snapshot] = fwdTimeCpmlFor2dAw(v, source, nDiffOrder, nBoundary, dz, dx, dt)
%
% FWDTIMECPMLFOR2DAW Simulate 2-d acoustic wave forward propagation using
% finite difference in time domain with the following partial differential
% equation (PDE)
%
% (1/v^2)*(d^2)u(z, x, t)/dt^2 + f(z, x, t) = zP + xP
% Update xPhi, zPhi, xA, zA, xPsi, zPsi, xP, zP and solve u(z, x, t) with Nonsplit Convolutional-PML (CPML)
%
% input arguments
% v(nz,nx)          velocity model
% source(nz,nx,nt)  source vector (e.g., shots)
% nDiffOrder        number of approximation order for differentiator operator
% nBoundary         thickness of the absorbing boundary
% dx                horizontal distance per sample
% dz                depth distance per sample
% dt                time difference per sample
%
% output arguments
% data              received 2-d x - time signal on the surface
% snapshot          pressure field u(z, x, t) in time domain
%
% Reference: 
% D. Komatitsch and R. Martin, An unsplit convolutional perfectly matched
% layer improved at grazing incidence for the seismic wave equation,
% Geophysics, Vol. 72 No. 5, pp. SM155-SM167, 2007
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com), Entao Liu (liuentao@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

%% Initialize storage
[nz, nx, nt] = size(source);
data = zeros(nx, nt);
coeff = dCoef(nDiffOrder, 's');
k = 2 * nDiffOrder - 1;
fdm = zeros(nz+2*k, nx+2*k, 3);  % fdm(:, :, 1) - past; fdm(:, :, 2) - present; fdm(:, :, 3) - future


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


%% 2-D Acoustic Wave Forward-Time Modeling

% finite difference coefficients
vdtSq = (v*dt).^2;

% common indicies
izi = (1+k):(nz+k);                         % used for z interior of fdm after zero-padding
ixi = (1+k):(nx+k);                         % used for x interior of fdm after zero-padding
izl = nDiffOrder:(nz+2*k-nDiffOrder);       % updated indices for zA
ixl = nDiffOrder:(nx+2*k-nDiffOrder);       % updated indices for xA

% snapshot
snapshot = zeros(nz, nx, nt);


% filenameFm2d = './videos/fm2d.avi';
% objVideoFm2d = VideoWriter(filenameFm2d);
% open(objVideoFm2d);

for it = 1:nt
    
    % finite differencing after zero-padding
    
    zPhi(izi, :) = zb .* zPhi(izi, :) + (zb - 1) .* diffOperator(fdm(izl+1, ixi, 2), coeff, dz, 1);
    xPhi(:, ixi) = xb .* xPhi(:, ixi) + (xb - 1) .* diffOperator(fdm(izi, ixl+1, 2), coeff, dx, 2);
    zA(izl, :) = diffOperator(fdm(:, ixi, 2), coeff, dz, 1) + zPhi(izl, :);
    xA(:, ixl) = diffOperator(fdm(izi, :, 2), coeff, dx, 2) + xPhi(:, ixl);
    
    zPsi(izi, :) = zb .* zPsi(izi, :) + (zb - 1) .* diffOperator(zA(izl, :), coeff, dz, 1);
    xPsi(:, ixi) = xb .* xPsi(:, ixi) + (xb - 1) .* diffOperator(xA(:, ixl), coeff, dx, 2);
    zP(izi, :) = diffOperator(zA(izl, :), coeff, dz, 1) + zPsi(izi, :);
    xP(:, ixi) = diffOperator(xA(:, ixl), coeff, dx, 2) + xPsi(:, ixi);
    
    fdm(izi, ixi, 3) = vdtSq .* (zP(izi, :) + xP(:, ixi) - source(:, :, it)) + 2 * fdm(izi, ixi, 2) - fdm(izi, ixi, 1);
    
    % *******************************************************************
    % FCT elimination of numerical dispersion. Needed for high frequency
    % source function and makes the algorithm very slow!!!)
    % *******************************************************************
    % fdm(:,:,3) = fctForAw(fdm(:,:,1),fdm(:,:,2),fdm(:,:,3));
    
    % update fdm for next time iteration
    fdm(:, :, 1) = fdm(:, :, 2);
    fdm(:, :, 2) = fdm(:, :, 3);
    
    % update data
    data(:, it) = fdm(1+k, ixi, 2);
    
    % update snapshot
    snapshot(:, :, it) = fdm(izi, ixi, 2);
    
%     % figure(3); imagesc(snapshot(1:nz-nBoundary,nBoundary+1:nx-nBoundary,it));
%     figure(3); imagesc(fdm(:, :, 2));
%     title(['Iteration: ',num2str(it)])
%     colorbar;
%     drawnow
%     writeVideo(objVideoFm2d, getframe(gcf));

end % time loop

% close(objVideoFm2d);

%data = data(21:nx-nBoundary,:);