function [data, snapshot] = fwdTimeCpmlFor3dAw(v, source, nDiffOrder, nBoundary, dz, dx, dy, dt)
%
% FWDTIMECPMLFOR3DAW Simulate 3-d acoustic wave forward propagating using
% finite difference in time domain with the following partial differential
% equation (PDE)
%
% (1/v^2)*(d^2)u(z, x, y, t)/dt^2 = zP + xP + yP + s(z, x, y, t)
% Update xPhi, yPhi, zPhi, xA, yA, zA, xPsi, yPsi, zPsi, xP, yP, zP and
% solve u(z, x, y, t) with Nonsplit Convolutional-PML (CPML)
%
% input arguments
% v(nz,nx,ny)           velocity model
% source(nz,nx,ny,nt)   source vector (e.g., shots)
% nDiffOrder            number of order of differentiator operator
% nBoundary             thickness of the absorbing boundary
% dx                    easting distance per sample
% dy                    northing distance per sample
% dz                    depth distance per sample
% dt                    time difference per sample
%
% output arguments
% data              received 3-d (x, y) - time signal on the surface
% snapshot          pressure field u(z, x, y, t) in time domain
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
[nz, nx, ny, nt] = size(source);
data = zeros(nx, ny, nt);
coeff = dCoef(nDiffOrder, 's');
k = 2 * nDiffOrder - 1;
fdm = zeros(nz+2*k, nx+2*k, ny+2*k, 3);  % fdm(:, :, :, 1) - past; fdm(:, :, :, 2) - present; fdm(:, :, :, 3) - future


%% Absorbing boundary condition (ABC): PML damping profile
ixb  = 1:nBoundary;         % index of x outside left boundary
ixb2 = nx-nBoundary+ixb;    % index of x outside right boundary
iyb  = 1:nBoundary;         % index of y outside front boundary
iyb2 = ny-nBoundary+iyb;    % index of y outside rear boundary
izb  = 1:nz-nBoundary;      % index of z inside down boundary
izb2 = nz-nBoundary+ixb;    % index of z outside down boundary

xDampLeft0 = dampPml(repmat(fliplr(ixb), nz, 1), v(:, ixb, 1), nBoundary);
xDampLeft = repmat(xDampLeft0, 1, 1, ny);
xDampRight0 = dampPml(repmat(ixb, nz, 1), v(:, ixb2, 1), nBoundary);
xDampRight = repmat(xDampRight0, 1, 1, ny);
xDamp = cat(2, xDampLeft, zeros(nz, nx-2*nBoundary, ny), xDampRight);
xb = exp(-xDamp * dt);

yDampFront0 = dampPml(repmat(reshape(fliplr(iyb), 1, 1, nBoundary), 1, nx, 1), v(1, :, iyb), nBoundary);
yDampFront = repmat(yDampFront0, nz, 1, 1);
yDampRear0 = dampPml(repmat(reshape(iyb, 1, 1, nBoundary), 1, nx, 1), v(1, :, iyb2), nBoundary);
yDampRear = repmat(yDampRear0, nz, 1, 1);
yDamp = cat(3, yDampFront, zeros(nz, nx, ny-2*nBoundary), yDampRear);
yb = exp(-yDamp * dt);

zDampDown0 = dampPml(repmat(ixb.', 1, nx), v(izb2, :, 1), nBoundary);
zDampDown = repmat(zDampDown0, 1, 1, ny);
zDamp = cat(1, zeros(nz-nBoundary, nx, ny), zDampDown);
zb = exp(-zDamp * dt);

% additional arrays for storage
zPhi = zeros(nz+2*k, nx, ny);
xPhi = zeros(nz, nx+2*k, ny);
yPhi = zeros(nz, nx, ny+2*k);
zA = zeros(nz+2*k, nx, ny);
xA = zeros(nz, nx+2*k, ny);
yA = zeros(nz, nx, ny+2*k);
zPsi = zeros(nz+k, nx, ny);
xPsi = zeros(nz, nx+k, ny);
yPsi = zeros(nz, nx, ny+k);
zP = zeros(nz+k, nx, ny);
xP = zeros(nz, nx+k, ny);
yP = zeros(nz, nx, ny+k);



%% 3-D Acoustic Wave Forward-Time Modeling

% finite difference coefficients
vdtSq = (v*dt).^2;

% common indicies
izi = (1+k):(nz+k);                         % used for z interior of fdm after zero-padding
ixi = (1+k):(nx+k);                         % used for x interior of fdm after zero-padding
iyi = (1+k):(ny+k);                         % used for y interior of fdm after zero-padding
izl = nDiffOrder:(nz+2*k-nDiffOrder);       % updated indices for zA
ixl = nDiffOrder:(nx+2*k-nDiffOrder);       % updated indices for xA
iyl = nDiffOrder:(ny+2*k-nDiffOrder);       % updated indices for xA


% snapshot
snapshot = zeros(nz, nx, ny, nt);


% filenameFm2d = './videos/fm2d.avi';
% objVideoFm2d = VideoWriter(filenameFm2d);
% open(objVideoFm2d);

for it = 1:nt
    
    % finite differencing after zero-padding
    
    zPhi(izi, :, :) = zb .* zPhi(izi, :, :) + (zb - 1) .* diffOperator(fdm(izl+1, ixi, iyi, 2), coeff, dz, 1);
    xPhi(:, ixi, :) = xb .* xPhi(:, ixi, :) + (xb - 1) .* diffOperator(fdm(izi, ixl+1, iyi, 2), coeff, dx, 2);
    yPhi(:, :, iyi) = yb .* yPhi(:, :, iyi) + (yb - 1) .* diffOperator(fdm(izi, ixi, iyl+1, 2), coeff, dy, 3);
    
    zA(izl, :, :) = diffOperator(fdm(:, ixi, iyi, 2), coeff, dz, 1) + zPhi(izl, :, :);
    xA(:, ixl, :) = diffOperator(fdm(izi, :, iyi, 2), coeff, dx, 2) + xPhi(:, ixl, :);
    yA(:, :, iyl) = diffOperator(fdm(izi, ixi, :, 2), coeff, dy, 3) + yPhi(:, :, iyl);
    
    zPsi(izi, :, :) = zb .* zPsi(izi, :, :) + (zb - 1) .* diffOperator(zA(izl, :, :), coeff, dz, 1);
    xPsi(:, ixi, :) = xb .* xPsi(:, ixi, :) + (xb - 1) .* diffOperator(xA(:, ixl, :), coeff, dx, 2);
    yPsi(:, :, iyi) = yb .* yPsi(:, :, iyi) + (yb - 1) .* diffOperator(yA(:, :, iyl), coeff, dy, 3);
    
    zP(izi, :, :) = diffOperator(zA(izl, :, :), coeff, dz, 1) + zPsi(izi, :, :);
    xP(:, ixi, :) = diffOperator(xA(:, ixl, :), coeff, dx, 2) + xPsi(:, ixi, :);
    yP(:, :, iyi) = diffOperator(yA(:, :, iyl), coeff, dy, 3) + yPsi(:, :, iyi);
    
    
    fdm(izi, ixi, iyi, 3) = vdtSq .* (zP(izi, :, :) + xP(:, ixi, :) + yP(:, :, iyi) + source(:, :, :, it)) + 2 * fdm(izi, ixi, iyi, 2) - fdm(izi, ixi, iyi, 1);
    
    % *******************************************************************
    % FCT elimination of numerical dispersion. Needed for high frequency
    % source function and makes the algorithm very slow!!!)
    % *******************************************************************
    % fdm(:,:,3) = fctForAw(fdm(:,:,1),fdm(:,:,2),fdm(:,:,3));
    
    % update fdm for next time iteration
    fdm(:, :, :, 1) = fdm(:, :, :, 2);
    fdm(:, :, :, 2) = fdm(:, :, :, 3);
    
    % update data
    data(:, :, it) = fdm(1+k, ixi, iyi, 2);
    
    % update snapshot
    snapshot(:, :, :, it) = fdm(izi, ixi, iyi, 2);
    
    
    
    %     % figure(3); imagesc(snapshot(1:nz-nBoundary,nBoundary+1:nx-nBoundary,it));
    %     figure(3); imagesc(fdm(:, :, 2));
    %     title(['Iteration: ',num2str(it)])
    %     colorbar;
    %     drawnow
    %     writeVideo(objVideoFm2d, getframe(gcf));
    
    
    
end % time loop

% close(objVideoFm2d);

%data = data(21:nx-nBoundary,:);