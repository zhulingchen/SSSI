function [snapshotVzp, snapshotVxp, snapshotVyp, snapshotVzs, snapshotVxs, snapshotVys] = fwdTimeSpmlFor3dEw(vp, vs, source, nDiffOrder, nBoundary, dz, dx, dy, dt)
%
% FWDTIMESPMLFOR3DEW Simulate 3-d elastic wave forward propagation using
% finite difference in time domain with the partial differential equations
% (PDEs) using split perfectly matched layer (SPML)
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

%% Initialize storage
[nz, nx, ny, nt] = size(source);

snapshotVzp = zeros(nz, nx, ny, nt);
snapshotVxp = zeros(nz, nx, ny, nt);
snapshotVyp = zeros(nz, nx, ny, nt);
snapshotVzs = zeros(nz, nx, ny, nt);
snapshotVxs = zeros(nz, nx, ny, nt);
snapshotVys = zeros(nz, nx, ny, nt);

coeff = dCoef(nDiffOrder, 's');
iz = 1+nDiffOrder:nz+nDiffOrder;      % interior z
ix = 1+nDiffOrder:nx+nDiffOrder;      % interior x
iy = 1+nDiffOrder:ny+nDiffOrder;      % interior y


%% Absorbing boundary condition (ABC): PML damping profile
ixb  = 1:nBoundary;         % index of x outside left boundary
ixb2 = nx-nBoundary+ixb;    % index of x outside right boundary
iyb  = 1:nBoundary;         % index of y outside front boundary
iyb2 = ny-nBoundary+iyb;    % index of y outside rear boundary
izb  = 1:nz-nBoundary;      % index of z inside down boundary
izb2 = nz-nBoundary+ixb;    % index of z outside down boundary

xDampLeft0 = dampPml(repmat(fliplr(ixb), nz, 1), vp(:, ixb, 1), nBoundary);
xDampLeft = repmat(xDampLeft0, 1, 1, ny);
xDampRight0 = dampPml(repmat(ixb, nz, 1), vp(:, ixb2, 1), nBoundary);
xDampRight = repmat(xDampRight0, 1, 1, ny);
xDamp = cat(2, xDampLeft, zeros(nz, nx-2*nBoundary, ny), xDampRight);
xb = exp(-xDamp * dt);

yDampFront0 = dampPml(repmat(reshape(fliplr(iyb), 1, 1, nBoundary), 1, nx, 1), vp(1, :, iyb), nBoundary);
yDampFront = repmat(yDampFront0, nz, 1, 1);
yDampRear0 = dampPml(repmat(reshape(iyb, 1, 1, nBoundary), 1, nx, 1), vp(1, :, iyb2), nBoundary);
yDampRear = repmat(yDampRear0, nz, 1, 1);
yDamp = cat(3, yDampFront, zeros(nz, nx, ny-2*nBoundary), yDampRear);
yb = exp(-yDamp * dt);

zDampDown0 = dampPml(repmat(ixb.', 1, nx), vp(izb2, :, 1), nBoundary);
zDampDown = repmat(zDampDown0, 1, 1, ny);
zDamp = cat(1, zeros(nz-nBoundary, nx, ny), zDampDown);
zb = exp(-zDamp * dt);

