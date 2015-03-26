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


%% Additional arrays for storage
% A is the divergence of the displacement vector S = (u, v, w).'
A_x = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
A_x = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
A_y = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
A_y = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
A_z = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
A_z = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
A = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
% B = (B1, B2, B3).' is the curl of the displacement vector S = (u, v, w).'
B1_y = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
B1_z = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
B1 = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
B2_z = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
B2_x = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
B2 = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
B3_x = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
B3_y = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
B3 = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);

Vxp = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
Vxs_y = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
Vxs_z = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
Vxs = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
Vx = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
Vyp = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
Vys_x = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
Vys_z = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
Vys = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
Vy = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
Vzp = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
Vzs_x = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
Vzs_y = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
Vzs = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);
Vz = zeros(nz+2*nDiffOrder, nx+2*nDiffOrder, ny+2*nDiffOrder, 3);


%% 3-D Elastic Wave Forward-Time Modeling

vpSq = vp.^2;
vpSq = padarray(vpSq, [nDiffOrder, nDiffOrder, nDiffOrder], 'replicate');
vsSq = vs.^2;
vsSq = padarray(vsSq, [nDiffOrder, nDiffOrder, nDiffOrder], 'replicate');

for it = 1:nt
    % calculate divergence and curl field of the displacement vector
    A_x(iz, ix, iy, 2) = (1 - 0.5 * xDamp * dt) ./ (1 + 0.5 * xDamp * dt) .* A_x(iz, ix, iy, 1) ...
        + dt ./ (1 + 0.5 * xDamp * dt) .* diffOperator(Vx(iz, 2:end, iy, 2), coeff, dx, 2);
    
    A_y(iz, ix, iy, 2) = (1 - 0.5 * yDamp * dt) ./ (1 + 0.5 * yDamp * dt) .* A_y(iz, ix, iy, 1) ...
        + dt ./ (1 + 0.5 * yDamp * dt) .* diffOperator(Vy(iz, ix, 1:end-1, 2), coeff, dy, 3);
    
    A_z(iz, ix, iy, 2) = (1 - 0.5 * zDamp * dt) ./ (1 + 0.5 * zDamp * dt) .* A_z(iz, ix, iy, 1) ...
        + dt ./ (1 + 0.5 * zDamp * dt) .* diffOperator(Vz(1:end-1, ix, iy, 2), coeff, dz, 1);
    
    A(iz, ix, iy, 2) = A_x(iz, ix, iy, 2) + A_y(iz, ix, iy, 2) + A_z(iz, ix, iy, 2);
    
    B1_y(iz, ix, iy, 2) = (1 - 0.5 * yDamp * dt) ./ (1 + 0.5 * yDamp * dt) .* B1_y(iz, ix, iy, 1) ...
        + dt ./ (1 + 0.5 * yDamp * dt) .* diffOperator(Vz(iz, ix, 2:end, 2), coeff, dy, 3);
    
    B1_z(iz, ix, iy, 2) = (1 - 0.5 * zDamp * dt) ./ (1 + 0.5 * zDamp * dt) .* B1_z(iz, ix, iy, 1) ...
        - dt ./ (1 + 0.5 * zDamp * dt) .* diffOperator(Vy(2:end, ix, iy, 2), coeff, dz, 1);
    
    B1(iz, ix, iy, 2) = B1_y(iz, ix, iy, 2) + B1_z(iz, ix, iy, 2);
    
    B2_z(iz, ix, iy, 2) = (1 - 0.5 * zDamp * dt) ./ (1 + 0.5 * zDamp * dt) .* B2_z(iz, ix, iy, 1) ...
        + dt ./ (1 + 0.5 * zDamp * dt) .* diffOperator(Vx(2:end, ix, iy, 2), coeff, dz, 1);
    
    B2_x(iz, ix, iy, 2) = (1 - 0.5 * xDamp * dt) ./ (1 + 0.5 * xDamp * dt) .* B2_x(iz, ix, iy, 1) ...
        - dt ./ (1 + 0.5 * xDamp * dt) .* diffOperator(Vz(iz, 1:end-1, iy, 2), coeff, dx, 2);
    
    B2(iz, ix, iy, 2) = B2_z(iz, ix, iy, 2) + B2_x(iz, ix, iy, 2);
    
    B3_x(iz, ix, iy, 2) = (1 - 0.5 * xDamp * dt) ./ (1 + 0.5 * xDamp * dt) .* B3_x(iz, ix, iy, 1) ...
        + dt ./ (1 + 0.5 * xDamp * dt) .* diffOperator(Vy(iz, 1:end-1, iy, 2), coeff, dx, 2);
    
    B3_y(iz, ix, iy, 2) = (1 - 0.5 * yDamp * dt) ./ (1 + 0.5 * yDamp * dt) .* B3_y(iz, ix, iy, 1) ...
        - dt ./ (1 + 0.5 * yDamp * dt) .* diffOperator(Vx(iz, ix, 2:end, 2), coeff, dy, 3);
    
    B3(iz, ix, iy, 2) = B3_x(iz, ix, iy, 2) + B3_y(iz, ix, iy, 2);
    
    % calculate velocity fields
    % P-wave field
    Vxp(iz, ix, iy, 3) = (1 - 0.5 * xDamp * dt) ./ (1 + 0.5 * xDamp * dt) .* Vxp(iz, ix, iy, 2) ...
        + vpSq(iz, ix, iy) * dt ./ (1 + 0.5 * xDamp * dt) .* diffOperator(A(iz, 1:end-1, iy, 2), coeff, dx, 2);
    
    Vyp(iz, ix, iy, 3) = (1 - 0.5 * yDamp * dt) ./ (1 + 0.5 * yDamp * dt) .* Vyp(iz, ix, iy, 2) ...
        + 0.25 * (vpSq(iz, ix, iy) + vpSq(iz, ix+1, iy) + vpSq(iz, ix+1, iy+1) + vpSq(iz, ix, iy+1)) * dt ./ (1 + 0.5 * yDamp * dt) .* diffOperator(A(iz, ix, 2:end, 2), coeff, dy, 3);
    
    Vzp(iz, ix, iy, 3) = (1 - 0.5 * zDamp * dt) ./ (1 + 0.5 * zDamp * dt) .* Vzp(iz, ix, iy, 2) ...
        + 0.25 * (vpSq(iz, ix, iy) + vpSq(iz, ix+1, iy) + vpSq(iz+1, ix+1, iy) + vpSq(iz+1, ix, iy)) * dt ./ (1 + 0.5 * zDamp * dt) .* diffOperator(A(2:end, ix, iy, 2), coeff, dz, 1);
    
    % S-wave field
    Vxs_y(iz, ix, iy, 3) = (1 - 0.5 * yDamp * dt) ./ (1 + 0.5 * yDamp * dt) .* Vxs_y(iz, ix, iy, 2) ...
        - 0.25 * (vsSq(iz-1, ix-1, iy) + vsSq(iz-1, ix, iy) + vsSq(iz, ix, iy) + vsSq(iz, ix-1, iy)) * dt ./ (1 + 0.5 * yDamp * dt) .* diffOperator(B3(iz, ix, 1:end-1, 2), coeff, dy, 3);
    
    Vxs_z(iz, ix, iy, 3) = (1 - 0.5 * zDamp * dt) ./ (1 + 0.5 * zDamp * dt) .* Vxs_z(iz, ix, iy, 2) ...
        + 0.25 * (vsSq(iz-1, ix-1, iy) + vsSq(iz-1, ix, iy) + vsSq(iz, ix, iy) + vsSq(iz, ix-1, iy)) * dt ./ (1 + 0.5 * zDamp * dt) .* diffOperator(B2(1:end-1, ix, iy, 2), coeff, dz, 1);
    
    Vxs(iz, ix, iy, 3) = Vxs_y(iz, ix, iy, 3) + Vxs_z(iz, ix, iy, 3) + source(:, :, :, it);      % source term
    
    Vys_x(iz, ix, iy, 3) = (1 - 0.5 * xDamp * dt) ./ (1 + 0.5 * xDamp * dt) .* Vys_x(iz, ix, iy, 2) ...
        + 0.25 * (vsSq(iz-1, ix, iy) + vsSq(iz-1, ix, iy+1) + vsSq(iz, ix, iy+1) + vsSq(iz, ix, iy)) * dt ./ (1 + 0.5 * xDamp * dt) .* diffOperator(B3(iz, 2:end, iy, 2), coeff, dx, 2);
    
    Vys_z(iz, ix, iy, 3) = (1 - 0.5 * zDamp * dt) ./ (1 + 0.5 * zDamp * dt) .* Vys_z(iz, ix, iy, 2) ...
        - 0.25 * (vsSq(iz-1, ix, iy) + vsSq(iz-1, ix, iy+1) + vsSq(iz, ix, iy+1) + vsSq(iz, ix, iy)) * dt ./ (1 + 0.5 * zDamp * dt) .* diffOperator(B1(1:end-1, ix, iy, 2), coeff, dz, 1);
    
    Vys(iz, ix, iy, 3) = Vys_x(iz, ix, iy, 3) + Vys_z(iz, ix, iy, 3) + source(:, :, :, it);      % source term
    
    Vzs_x(iz, ix, iy, 3) = (1 - 0.5 * xDamp * dt) ./ (1 + 0.5 * xDamp * dt) .* Vzs_x(iz, ix, iy, 2) ...
        - vsSq(iz, ix, iy) * dt ./ (1 + 0.5 * xDamp * dt) .* diffOperator(B2(iz, 2:end, iy, 2), coeff, dx, 2);
    
    Vzs_y(iz, ix, iy, 3) = (1 - 0.5 * yDamp * dt) ./ (1 + 0.5 * yDamp * dt) .* Vzs_y(iz, ix, iy, 2) ...
        + vsSq(iz, ix, iy) * dt ./ (1 + 0.5 * yDamp * dt) .* diffOperator(B1(iz, ix, 1:end-1, 2), coeff, dy, 3);
    
    Vzs(iz, ix, iy, 3) = Vzs_x(iz, ix, iy, 3) + Vzs_y(iz, ix, iy, 3) + source(:, :, :, it);      % source term
    
    %% velocity wavefield of each direction is composed by its P-wave and S-wave components
    Vx(iz, ix, iy, 3) = Vxp(iz, ix, iy, 3) + Vxs(iz, ix, iy, 3);
    Vy(iz, ix, iy, 3) = Vyp(iz, ix, iy, 3) + Vys(iz, ix, iy, 3);
    Vz(iz, ix, iy, 3) = Vzp(iz, ix, iy, 3) + Vzs(iz, ix, iy, 3);
    
    %% update the wavefield
    A_z(iz, ix, iy, 1) = A_z(iz, ix, iy, 2); A_z(iz, ix, iy, 2) = A_z(iz, ix, iy, 3);
    A_x(iz, ix, iy, 1) = A_x(iz, ix, iy, 2); A_x(iz, ix, iy, 2) = A_x(iz, ix, iy, 3);
    A(iz, ix, iy, 1) = A(iz, ix, iy, 2); A(iz, ix, iy, 2) = A(iz, ix, iy, 3);
    B1_y(iz, ix, iy, 1) = B1_y(iz, ix, iy, 2); B1_y(iz, ix, iy, 2) = B1_y(iz, ix, iy, 3);
    B1_z(iz, ix, iy, 1) = B1_z(iz, ix, iy, 2); B1_z(iz, ix, iy, 2) = B1_z(iz, ix, iy, 3);
    B1(iz, ix, iy, 1) = B1(iz, ix, iy, 2); B1(iz, ix, iy, 2) = B1(iz, ix, iy, 3);
    B2_z(iz, ix, iy, 1) = B2_z(iz, ix, iy, 2); B2_z(iz, ix, iy, 2) = B2_z(iz, ix, iy, 3);
    B2_x(iz, ix, iy, 1) = B2_x(iz, ix, iy, 2); B2_x(iz, ix, iy, 2) = B2_x(iz, ix, iy, 3);
    B2(iz, ix, iy, 1) = B2(iz, ix, iy, 2); B2(iz, ix, iy, 2) = B2(iz, ix, iy, 3);
    B3_x(iz, ix, iy, 1) = B3_x(iz, ix, iy, 2); B3_x(iz, ix, iy, 2) = B3_x(iz, ix, iy, 3);
    B3_y(iz, ix, iy, 1) = B3_y(iz, ix, iy, 2); B3_y(iz, ix, iy, 2) = B3_y(iz, ix, iy, 3);
    B3(iz, ix, iy, 1) = B3(iz, ix, iy, 2); B3(iz, ix, iy, 2) = B3(iz, ix, iy, 3);
    
    Vxp(iz, ix, iy, 1) = Vxp(iz, ix, iy, 2); Vxp(iz, ix, iy, 2) = Vxp(iz, ix, iy, 3);
    Vyp(iz, ix, iy, 1) = Vyp(iz, ix, iy, 2); Vyp(iz, ix, iy, 2) = Vyp(iz, ix, iy, 3);
    Vzp(iz, ix, iy, 1) = Vzp(iz, ix, iy, 2); Vzp(iz, ix, iy, 2) = Vzp(iz, ix, iy, 3);
    Vxs_y(iz, ix, iy, 1) = Vxs_y(iz, ix, iy, 2); Vxs_y(iz, ix, iy, 2) = Vxs_y(iz, ix, iy, 3);
    Vxs_z(iz, ix, iy, 1) = Vxs_z(iz, ix, iy, 2); Vxs_z(iz, ix, iy, 2) = Vxs_z(iz, ix, iy, 3);
    Vxs(iz, ix, iy, 1) = Vxs(iz, ix, iy, 2); Vxs(iz, ix, iy, 2) = Vxs(iz, ix, iy, 3);
    Vys_x(iz, ix, iy, 1) = Vys_x(iz, ix, iy, 2); Vys_x(iz, ix, iy, 2) = Vys_x(iz, ix, iy, 3);
    Vys_z(iz, ix, iy, 1) = Vys_z(iz, ix, iy, 2); Vys_z(iz, ix, iy, 2) = Vys_z(iz, ix, iy, 3);
    Vys(iz, ix, iy, 1) = Vys(iz, ix, iy, 2); Vys(iz, ix, iy, 2) = Vys(iz, ix, iy, 3);
    Vzs_x(iz, ix, iy, 1) = Vzs_x(iz, ix, iy, 2); Vzs_x(iz, ix, iy, 2) = Vzs_x(iz, ix, iy, 3);
    Vzs_y(iz, ix, iy, 1) = Vzs_y(iz, ix, iy, 2); Vzs_y(iz, ix, iy, 2) = Vzs_y(iz, ix, iy, 3);
    Vzs(iz, ix, iy, 1) = Vzs(iz, ix, iy, 2); Vzs(iz, ix, iy, 2) = Vzs(iz, ix, iy, 3);
    
    Vx(iz, ix, iy, 1) = Vx(iz, ix, iy, 2); Vx(iz, ix, iy, 2) = Vx(iz, ix, iy, 3);
    Vy(iz, ix, iy, 1) = Vy(iz, ix, iy, 2); Vy(iz, ix, iy, 2) = Vy(iz, ix, iy, 3);
    Vz(iz, ix, iy, 1) = Vz(iz, ix, iy, 2); Vz(iz, ix, iy, 2) = Vz(iz, ix, iy, 3);
    
    %% P-wave and S-wave snapshots of x-axis and z-axis velocity wavefield
    snapshotVxp(:, :, :, it) = Vxp(iz, ix, iy, 2);
    snapshotVyp(:, :, :, it) = Vyp(iz, ix, iy, 2);
    snapshotVzp(:, :, :, it) = Vzp(iz, ix, iy, 2);
    snapshotVxs(:, :, :, it) = Vxs(iz, ix, iy, 2);
    snapshotVys(:, :, :, it) = Vys(iz, ix, iy, 2);
    snapshotVzs(:, :, :, it) = Vzs(iz, ix, iy, 2);
    
end



