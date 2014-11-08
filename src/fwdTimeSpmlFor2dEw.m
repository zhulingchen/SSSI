function [dataVzp, dataVxp, dataVzs, dataVxs] = fwdTimeSpmlFor2dEw(vp, vs, source, nDiffOrder, nBoundary, dz, dx, dt)
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

dataVzp = zeros(nz, nx, nt);
dataVxp = zeros(nz, nx, nt);
dataVzs = zeros(nz, nx, nt);
dataVxs = zeros(nz, nx, nt);

coeff = dCoef(nDiffOrder, 's');
iz = 1+nDiffOrder:nz-nDiffOrder;      % interior z
ix = 1+nDiffOrder:nx-nDiffOrder;      % interior x

%% Absorbing boundary condition (ABC): Split PML (SPML)
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
A_perp = zeros(nz, nx, 3);
A_para = zeros(nz, nx, 3);
B_perp = zeros(nz, nx, 3);
B_para = zeros(nz, nx, 3);
A = zeros(nz, nx, 3);
B = zeros(nz, nx, 3);

Vxp = zeros(nz, nx, 3);
Vxs = zeros(nz, nx, 3);
Vx = zeros(nz, nx, 3);
Vzp = zeros(nz, nx, 3);
Vzs = zeros(nz, nx, 3);
Vz = zeros(nz, nx, 3);

dataVx = zeros(nz, nx, nt);
dataVz = zeros(nz, nx, nt);

%% 2-D Elastic Wave Forward-Time Modeling

vpSq = vp.^2;
vsSq = vs.^2;

for it = 1:nt
    A_perp(iz, ix, 2) = (1 - 0.5 * zDamp(iz, ix) * dt) ./ (1 + 0.5 * zDamp(iz, ix) * dt) .* A_perp(iz, ix, 1) ...
        + dt ./ (1 + 0.5 * zDamp(iz, ix) * dt) .* diffOperator(Vz(1:end-1, ix, 2), coeff, dz, 1);
    
    A_para(iz, ix, 2) = (1 - 0.5 * xDamp(iz, ix) * dt) ./ (1 + 0.5 * xDamp(iz, ix) * dt) .* A_para(iz, ix, 1) ...
        + dt ./ (1 + 0.5 * xDamp(iz, ix) * dt) .* diffOperator(Vx(iz, 2:end, 2), coeff, dx, 2);
    
    B_perp(iz, ix, 2) = (1 - 0.5 * zDamp(iz, ix) * dt) ./ (1 + 0.5 * zDamp(iz, ix) * dt) .* B_perp(iz, ix, 1)...
        + dt ./ (1 + 0.5 * zDamp(iz, ix) * dt) .* diffOperator(Vx(2:end, ix, 2), coeff, dz, 1);
    
    B_para(iz, ix, 2) = (1 - 0.5 * xDamp(iz, ix) * dt) ./ (1 + 0.5 * xDamp(iz, ix) * dt) .* B_para(iz, ix, 1)...
        - dt ./ (1 + 0.5 * xDamp(iz, ix) * dt) .* diffOperator(Vz(iz, 1:end-1, 2), coeff, dx, 2);
    
    A(iz, ix, 2) = A_perp(iz, ix, 2) + A_para(iz, ix, 2);
    
    B(iz, ix, 2) = B_perp(iz, ix, 2) + B_para(iz, ix, 2);
    
    Vzp(iz, ix, 3) = (1 - 0.5 * zDamp(iz, ix) * dt) ./ (1 + 0.5 * zDamp(iz, ix) * dt) .* Vzp(iz, ix, 2) ...
        + 0.25 .* (vpSq(iz+1, ix) + vpSq(iz+1, ix+1) + vpSq(iz, ix) + vpSq(iz, ix+1)) * dt ./ (1 + 0.5 * zDamp(iz, ix) * dt) ...
        .* diffOperator(A(2:end, ix, 2), coeff, dz, 1) ...
        + source(iz, ix, it);      % soure term
    
    Vxp(iz, ix, 3) = (1 - 0.5 * xDamp(iz, ix) * dt) ./ (1 + 0.5 * xDamp(iz, ix) * dt) .* Vxp(iz, ix, 2) ...
        + vpSq(iz, ix) * dt ./ (1 + 0.5 * xDamp(iz, ix) * dt) ...
        .* diffOperator(A(iz, 1:end-1, 2), coeff, dx, 2);
    
    Vzs(iz, ix, 3) = (1 - 0.5 * xDamp(iz, ix) * dt) ./ (1 + 0.5 * xDamp(iz, ix) * dt) .* Vzs(iz, ix, 2) ...
        - 0.25 .* (vsSq(iz, ix) + vsSq(iz, ix-1) + vsSq(iz-1, ix) + vsSq(iz-1, ix-1)) * dt ./ (1 + 0.5 * xDamp(iz, ix) * dt) ...
        .* diffOperator(B(iz, 2:end, 2), coeff, dx, 2) ...
        + source(iz, ix, it);      % source term
    
    Vxs(iz, ix, 3) = (1 - 0.5 * zDamp(iz, ix) * dt) ./ (1 + 0.5 * zDamp(iz, ix) * dt) .* Vxs(iz, ix, 2)...
        + vsSq(iz, ix) * dt ./ (1 + 0.5 * zDamp(iz, ix) * dt) ...
        .* diffOperator(B(1:end-1, ix, 2), coeff, dz, 1);
    
    
    %% velocity wavefield of each direction is composed by its P-wave and S-wave components
    Vx(iz, ix, 3) = Vxp(iz, ix, 3) + Vxs(iz, ix, 3);
    
    Vz(iz, ix, 3) = Vzp(iz, ix, 3) + Vzs(iz, ix, 3);
    
    %% record snapshot for x-axis component and z-axis component
    dataVx(:, :, it) = Vx(:, :, 3);
    dataVz(:, :, it) = Vz(:, :, 3);
    dataVxp(:, :, it) = Vxp(:, :, 3);
    dataVzp(:, :, it) = Vzp(:, :, 3);
    dataVxs(:, :, it) = Vxs(:, :, 3);
    dataVzs(:, :, it) = Vzs(:, :, 3);
    
    %% update the wavefield
    A_perp(iz, ix, 1) = A_perp(iz, ix, 2); A_perp(iz, ix, 2) = A_perp(iz, ix, 3);
    A_para(iz, ix, 1) = A_para(iz, ix, 2); A_para(iz, ix, 2) = A_para(iz, ix, 3);
    B_perp(iz, ix, 1) = B_perp(iz, ix, 2); B_perp(iz, ix, 2) = B_perp(iz, ix, 3);
    B_para(iz, ix, 1) = B_para(iz, ix, 2); B_para(iz, ix, 2) = B_para(iz, ix, 3);
    
    A(iz, ix, 1) = A(iz, ix, 2); A(iz, ix, 2) = A(iz, ix, 3);
    B(iz, ix, 1) = B(iz, ix, 2); B(iz, ix, 2) = B(iz, ix, 3);
    
    Vxp(iz, ix, 1) = Vxp(iz, ix, 2); Vxp(iz, ix, 2) = Vxp(iz, ix, 3);
    Vzp(iz, ix, 1) = Vzp(iz, ix, 2); Vzp(iz, ix, 2) = Vzp(iz, ix, 3);
    Vxs(iz, ix, 1) = Vxs(iz, ix, 2); Vxs(iz, ix, 2) = Vxs(iz, ix, 3);
    Vzs(iz, ix, 1) = Vzs(iz, ix, 2); Vzs(iz, ix, 2) = Vzs(iz, ix, 3);
    
    Vx(iz, ix, 1) = Vx(iz, ix, 2); Vx(iz, ix, 2) = Vx(iz, ix, 3);
    Vz(iz, ix, 1) = Vz(iz, ix, 2); Vz(iz, ix, 2) = Vz(iz, ix, 3);
    
end  % time loop ends
