function [snapshotVzp, snapshotVxp, snapshotVzs, snapshotVxs] = fwdTimeSpmlFor2dEw(vp, vs, source, nDiffOrder, nBoundary, dz, dx, dt)
%
% FWDTIMESPMLFOR2DEW Simulate 2-d elastic wave forward propagation using
% finite difference in time domain with the partial differential equations
% (PDEs) using split perfectly matched layer (SPML) from the following
% reference
%
% Ke-Yang Chen, Finite-Difference Simulation of Elastic Wave with
% Separation in Pure P- and S-Modes, Journal of Computational Methods in
% Physics, vol. 2014, Article ID 108713, 14 pages, 2014.
% doi:10.1155/2014/108713
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com), Entao Liu (liuentao@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

%% Initialize storage
[nz, nx, nt] = size(source);

snapshotVzp = zeros(nz, nx, nt);
snapshotVxp = zeros(nz, nx, nt);
snapshotVzs = zeros(nz, nx, nt);
snapshotVxs = zeros(nz, nx, nt);

coeff = dCoef(nDiffOrder, 's');
iz = 1+nDiffOrder:nz+nDiffOrder;      % interior z
ix = 1+nDiffOrder:nx+nDiffOrder;      % interior x

%% Absorbing boundary condition (ABC): PML damping profile
ixb = 1:nBoundary;          % index of x outside left boundary
ixb2 = nx-nBoundary+ixb;    % index of x outside right boundary
izb  = 1:nz-nBoundary;      % index of z inside down boundary
izb2 = nz-nBoundary+ixb;    % index of z outside down boundary

xDampLeft = dampPml(repmat(fliplr(ixb), nz, 1), vp(:, ixb), nBoundary);
xDampRight = dampPml(repmat(ixb, nz, 1), vp(:, ixb2), nBoundary);
xDamp = [xDampLeft, zeros(nz, nx-2*nBoundary), xDampRight];

zDampDown = dampPml(repmat(ixb.', 1, nx), vp(izb2, :), nBoundary);
zDamp = [zeros(nz-nBoundary, nx); zDampDown];

%% additional arrays for storage
% A is the divergence of the displacement vector S = (w, u).'
A_z = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, 3);
A_x = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, 3);
A = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, 3);
% B is the curl of the displacement vector S = (w, u).', whose direction
% points to the y-axis
B_z = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, 3);
B_x = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, 3);
B = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, 3);

Vxp = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, 3);
Vxs = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, 3);
Vx = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, 3);
Vzp = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, 3);
Vzs = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, 3);
Vz = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, 3);


%% 2-D Elastic Wave Forward-Time Modeling

vpSq = vp.^2;
vpSq = padarray(vpSq, [nDiffOrder, nDiffOrder], 'replicate');
vsSq = vs.^2;
vsSq = padarray(vsSq, [nDiffOrder, nDiffOrder], 'replicate');

for it = 1:nt
    
    A_x(iz, ix, 2) = (1 - 0.5 * xDamp * dt) ./ (1 + 0.5 * xDamp * dt) .* A_x(iz, ix, 1) ...
        + dt ./ (1 + 0.5 * xDamp * dt) .* diffOperator(Vx(iz, 1:end-1, 2), coeff, dx, 2);
    
    A_z(iz, ix, 2) = (1 - 0.5 * zDamp * dt) ./ (1 + 0.5 * zDamp * dt) .* A_z(iz, ix, 1) ...
        + dt ./ (1 + 0.5 * zDamp * dt) .* diffOperator(Vz(1:end-1, ix, 2), coeff, dz, 1);
    
    A(iz, ix, 2) = A_x(iz, ix, 2) + A_z(iz, ix, 2);
    
    B_z(iz, ix, 2) = (1 - 0.5 * zDamp * dt) ./ (1 + 0.5 * zDamp * dt) .* B_z(iz, ix, 1)...
        + dt ./ (1 + 0.5 * zDamp * dt) .* diffOperator(Vx(2:end, ix, 2), coeff, dz, 1);
    
    B_x(iz, ix, 2) = (1 - 0.5 * xDamp * dt) ./ (1 + 0.5 * xDamp * dt) .* B_x(iz, ix, 1)...
        - dt ./ (1 + 0.5 * xDamp * dt) .* diffOperator(Vz(iz, 2:end, 2), coeff, dx, 2);
    
    B(iz, ix, 2) = B_z(iz, ix, 2) + B_x(iz, ix, 2);
    
    Vzp(iz, ix, 3) = (1 - 0.5 * zDamp * dt) ./ (1 + 0.5 * zDamp * dt) .* Vzp(iz, ix, 2) ...
        + 0.25 .* (vpSq(iz, ix) + vpSq(iz+1, ix) + vpSq(iz+1, ix-1) + vpSq(iz, ix-1)) * dt ./ (1 + 0.5 * zDamp * dt) .* diffOperator(A(2:end, ix, 2), coeff, dz, 1) ...
        + dt ./ (1 + 0.5 * zDamp * dt) .* source(:, :, it);      % source term
    
    Vxp(iz, ix, 3) = (1 - 0.5 * xDamp * dt) ./ (1 + 0.5 * xDamp * dt) .* Vxp(iz, ix, 2) ...
        + vpSq(iz, ix) * dt ./ (1 + 0.5 * xDamp * dt) .* diffOperator(A(iz, 2:end, 2), coeff, dx, 2) ...
        + dt ./ (1 + 0.5 * xDamp * dt) .* source(:, :, it);      % source term
    
    Vzs(iz, ix, 3) = (1 - 0.5 * xDamp * dt) ./ (1 + 0.5 * xDamp * dt) .* Vzs(iz, ix, 2) ...
        - vsSq(iz, ix) * dt ./ (1 + 0.5 * xDamp * dt) ...
        .* diffOperator(B(iz, 1:end-1, 2), coeff, dx, 2);
    
    Vxs(iz, ix, 3) = (1 - 0.5 * zDamp * dt) ./ (1 + 0.5 * zDamp * dt) .* Vxs(iz, ix, 2)...
        + 0.25 .* (vsSq(iz, ix) + vsSq(iz, ix+1) + vsSq(iz-1, ix+1) + vsSq(iz-1, ix)) * dt ./ (1 + 0.5 * zDamp * dt) ...
        .* diffOperator(B(1:end-1, ix, 2), coeff, dz, 1);
    
    
    %% velocity wavefield of each direction is composed by its P-wave and S-wave components
    Vx(iz, ix, 3) = Vxp(iz, ix, 3) + Vxs(iz, ix, 3);
    Vz(iz, ix, 3) = Vzp(iz, ix, 3) + Vzs(iz, ix, 3);
    
    %% update the wavefield
    A_z(:, :, 1) = A_z(:, :, 2); A_z(:, :, 2) = A_z(:, :, 3);
    A_x(:, :, 1) = A_x(:, :, 2); A_x(:, :, 2) = A_x(:, :, 3);
    A(:, :, 1) = A(:, :, 2); A(:, :, 2) = A(:, :, 3);
    B_z(:, :, 1) = B_z(:, :, 2); B_z(:, :, 2) = B_z(:, :, 3);
    B_x(:, :, 1) = B_x(:, :, 2); B_x(:, :, 2) = B_x(:, :, 3);
    B(:, :, 1) = B(:, :, 2); B(:, :, 2) = B(:, :, 3);
    
    Vxp(:, :, 1) = Vxp(:, :, 2); Vxp(:, :, 2) = Vxp(:, :, 3);
    Vzp(:, :, 1) = Vzp(:, :, 2); Vzp(:, :, 2) = Vzp(:, :, 3);
    Vxs(:, :, 1) = Vxs(:, :, 2); Vxs(:, :, 2) = Vxs(:, :, 3);
    Vzs(:, :, 1) = Vzs(:, :, 2); Vzs(:, :, 2) = Vzs(:, :, 3);
    
    Vx(:, :, 1) = Vx(:, :, 2); Vx(:, :, 2) = Vx(:, :, 3);
    Vz(:, :, 1) = Vz(:, :, 2); Vz(:, :, 2) = Vz(:, :, 3);
    
    %% P-wave and S-wave snapshots of x-axis and z-axis velocity wavefield
    snapshotVxp(:, :, it) = Vxp(iz, ix, 2);
    snapshotVzp(:, :, it) = Vzp(iz, ix, 2);
    snapshotVxs(:, :, it) = Vxs(iz, ix, 2);
    snapshotVzs(:, :, it) = Vzs(iz, ix, 2);
    
end  % time loop ends
