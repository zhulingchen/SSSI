function [dataVxp, dataVzp, dataVxs, dataVzs] = fwdTimeSpmlFor2dEw(vp, vs, source, nDiffOrder, nBoundary, dz, dx, dt)
%
% FWDTIMESPMLFOR2DEW Simulate 2-d elastic wave forward propagation using
% finite difference in time domain with the partial differential equations
% (PDEs) using split perfectly matched layer (SPML) from the following
% reference
%
% Ke-Yang Chen, ?Finite-Difference Simulation of Elastic Wave with
% Separation in Pure P- and S-Modes,? Journal of Computational Methods in
% Physics, vol. 2014, Article ID 108713, 14 pages, 2014.
% doi:10.1155/2014/108713
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

%% Initialize storage
[nz, nx, nt] = size(source);
dataVxp = zeros(nz,nx,nt);
dataVzp = zeros(nz,nx,nt);
dataVxs = zeros(nz,nx,nt);
dataVzs = zeros(nz,nx,nt);

coeff = dCoef(nDiffOrder, 's');
iz = 1+nDiffOrder:nz-nDiffOrder;      % interior z
ix = 1+nDiffOrder:nx-nDiffOrder;      % interior x

% % common indicies
% k = 2 * nDiffOrder - 1;
% izi = (1+k):(nz+k);                         % used for z interior of fdm after zero-padding
% ixi = (1+k):(nx+k);                         % used for x interior of fdm after zero-padding
% izl = nDiffOrder:(nz+2*k-nDiffOrder);       % updated indices for zA
% ixl = nDiffOrder:(nx+2*k-nDiffOrder);       % updated indices for xA

%% Absorbing boundary condition (ABC): Split PML (SPML)
ixb = 1:nBoundary;          % index of x outside left boundary
ixb2 = nx-nBoundary+ixb;    % index of x outside right boundary
izb  = 1:nz-nBoundary;      % index of z inside down boundary
izb2 = nz-nBoundary+ixb;    % index of z outside down boundary

xDampLeft = dampPml(repmat(fliplr(ixb), nz, 1), vp(:, ixb), nBoundary);
xDampRight = dampPml(repmat(ixb, nz, 1), vp(:, ixb2), nBoundary);
xDamp = [xDampLeft, zeros(nz, nx-2*nBoundary), xDampRight];

zDampUp = dampPml(repmat(fliplr(ixb.'), 1, nx), vp(1:nBoundary, :), nBoundary);
zDampDown = dampPml(repmat(ixb.', 1, nx), vp(izb2, :), nBoundary);
zDamp = [zeros(nz-nBoundary, nx); zDampDown];

%% additional arrays for storage
Vx  = zeros(nz,nx,3);
Vxp = zeros(nz,nx,3);
Vxs = zeros(nz,nx,3);
Vz  = zeros(nz,nx,3);
Vzp = zeros(nz,nx,3);
Vzs = zeros(nz,nx,3);
A   = zeros(nz,nx,3);
B   = zeros(nz,nx,3);
A_1   = zeros(nz,nx,3);
A_2   = zeros(nz,nx,3);
B_1   = zeros(nz,nx,3);
B_2   = zeros(nz,nx,3);
dataVx = zeros(nz,nx,nt);
dataVz = zeros(nz,nx,nt);

%% Forward-Time Modeling
vp2 = vp.^2;
vs2 = vs.^2;
TOX = dt / dx;          % dt over dx
TOZ = dt / dz;          % dt over dz

for it=1:nt
    % In the reference the first index is referred to z direction and
    % second index is referred to x direction. So swap the xDamp and zDamp
    A_1(iz,ix,2) = (1-0.5.*zDamp(iz,ix).*dt)./(1+0.5.*zDamp(iz,ix).*dt).*A_1(iz,ix,1)...
        + TOX./(1+0.5.*zDamp(iz,ix).*dt) .* ...
        (coeff(1) * (Vx(iz+1,ix,2)-Vx(iz,ix,2)) ...
        + coeff(2) * (Vx(iz+2,ix,2)-Vx(iz-1,ix,2)) ...
        + coeff(3) * (Vx(iz+3,ix,2)-Vx(iz-2,ix,2)));
    
%     temp1 = (coeff(1) * (Vx(iz+1,ix,2)-Vx(iz,ix,2)) ...
%         + coeff(2) * (Vx(iz+2,ix,2)-Vx(iz-1,ix,2)) ...
%         + coeff(3) * (Vx(iz+3,ix,2)-Vx(iz-2,ix,2)));
%     
%     temp2 = diffOperator(Vx(iz+1, ix, 2), coeff, dz, 2) * dz;
    
    
    A_2(iz,ix,2) = (1-0.5.*xDamp(iz,ix).*dt)./(1+0.5.*xDamp(iz,ix).*dt).*A_2(iz,ix,1)...
        + TOZ./(1+0.5.*xDamp(iz,ix).*dt) .* ...
        (coeff(1) * (Vz(iz,ix,2)-Vz(iz,ix-1,2)) ...
        + coeff(2) * (Vz(iz,ix+1,2)-Vz(iz,ix-2,2)) ...
        + coeff(3) * (Vz(iz,ix+2,2)-Vz(iz,ix-3,2)));
    
    B_1(iz,ix,2) = (1-0.5.*zDamp(iz,ix).*dt)./(1+0.5.*zDamp(iz,ix).*dt).*B_1(iz,ix,1)...
        - TOX./(1+0.5.*zDamp(iz,ix).*dt) .* ...
        (coeff(1) * (Vz(iz,ix,2)-Vz(iz-1,ix,2)) ...
        + coeff(2) * (Vz(iz+1,ix,2)-Vz(iz-2,ix,2)) ...
        + coeff(3) * (Vz(iz+2,ix,2)-Vz(iz-3,ix,2)));
    
    B_2(iz,ix,2) = (1-0.5.*xDamp(iz,ix).*dt)./(1+0.5.*xDamp(iz,ix).*dt).*B_2(iz,ix,1)...
        + TOZ./(1+0.5.*xDamp(iz,ix).*dt) .* ...
        (coeff(1) * (Vx(iz,ix+1,2)-Vx(iz,ix,2)) ...
        + coeff(2) * (Vx(iz,ix+2,2)-Vx(iz,ix-1,2)) ...
        + coeff(3) * (Vx(iz,ix+3,2)-Vx(iz,ix-2,2)));
    
    A(iz,ix,2) = A_1(iz,ix,2) + A_2(iz,ix,2);
    
    B(iz,ix,2) = B_1(iz,ix,2) + B_2(iz,ix,2);
    
    Vxp(iz,ix,3) = (1-0.5.*zDamp(iz,ix).*dt)./(1+0.5.*zDamp(iz,ix).*dt).*Vxp(iz,ix,2)...
        + vp2(iz,ix).*TOX./(1+0.5.*zDamp(iz,ix).*dt) .* ...
        (coeff(1) * (A(iz,ix,2)-A(iz-1,ix,2)) ...
        + coeff(2) * (A(iz+1,ix,2)-A(iz-2,ix,2)) ...
        + coeff(3) * (A(iz+2,ix,2)-A(iz-3,ix,2))) ...
        + source(iz,ix,it);      % soure term
    
    Vzp(iz,ix,3) = (1-0.5.*xDamp(iz,ix).*dt)./(1+0.5.*xDamp(iz,ix).*dt).*Vzp(iz,ix,2)...
        + 0.25.*(vp2(iz+1,ix)+vp2(iz+1,ix+1)+vp2(iz,ix)+vp2(iz,ix+1))...
        .*TOZ./(1+0.5.*xDamp(iz,ix).*dt) .* ...
        (coeff(1) * (A(iz,ix+1,2)-A(iz,ix,2)) ...
        + coeff(2) * (A(iz,ix+2,2)-A(iz,ix-1,2)) ...
        + coeff(3) * (A(iz,ix+3,2)-A(iz,ix-2,2)));
    
    Vxs(iz,ix,3) = (1-0.5.*xDamp(iz,ix).*dt)./(1+0.5.*xDamp(iz,ix).*dt).*Vxs(iz,ix,2)...
        + vs2(iz,ix).*TOZ./(1+0.5.*xDamp(iz,ix).*dt) .* ...
        (coeff(1) * (B(iz,ix,2)-B(iz,ix-1,2)) ...
        + coeff(2) * (B(iz,ix+1,2)-B(iz,ix-2,2)) ...
        + coeff(3) * (B(iz,ix+2,2)-B(iz,ix-3,2))) ...
        + source(iz,ix,it);      % source term
    
    
    Vzs(iz,ix,3) = (1-0.5.*zDamp(iz,ix).*dt)./(1+0.5.*zDamp(iz,ix).*dt).*Vzs(iz,ix,2)...
        - 0.25.*(vs2(iz,ix)+vs2(iz,ix-1)+vs2(iz-1,ix)+vs2(iz-1,ix-1))...
        .*TOX./(1+0.5.*zDamp(iz,ix).*dt) .* ...
        (coeff(1) * (B(iz+1,ix,2)-B(iz,ix,2)) ...
        + coeff(2) * (B(iz+2,ix,2)-B(iz-1,ix,2)) ...
        + coeff(3) * (B(iz+3,ix,2)-B(iz-2,ix,2)));
    
    % The velocity wavefield of each direction is composed by its x-axis and z-axis components
    
    Vx(iz,ix,3) = Vxp(iz,ix,3) + Vxs(iz,ix,3);
    
    Vz(iz,ix,3) = Vzp(iz,ix,3) + Vzs(iz,ix,3);
    
    % ***********  record snapshot for x-axis component and z-axis component **************************
    dataVx(:,:,it) = Vx(:,:,3);
    dataVz(:,:,it) = Vz(:,:,3);
    dataVxp(:,:,it) = Vxp(:,:,3);
    dataVzp(:,:,it) = Vzp(:,:,3);
    dataVxs(:,:,it) = Vxs(:,:,3);
    dataVzs(:,:,it) = Vzs(:,:,3);
    
    % ***********  update the wavefield ******************************
    Vx(iz,ix,1) = Vx(iz,ix,2); Vx(iz,ix,2) = Vx(iz,ix,3);
    Vz(iz,ix,1) = Vz(iz,ix,2); Vz(iz,ix,2) = Vz(iz,ix,3);
    Vxp(iz,ix,1) = Vxp(iz,ix,2); Vxp(iz,ix,2) = Vxp(iz,ix,3);
    Vzp(iz,ix,1) = Vzp(iz,ix,2); Vzp(iz,ix,2) = Vzp(iz,ix,3);
    Vxs(iz,ix,1) = Vxs(iz,ix,2); Vxs(iz,ix,2) = Vxs(iz,ix,3);
    Vzs(iz,ix,1) = Vzs(iz,ix,2); Vzs(iz,ix,2) = Vzs(iz,ix,3);
    A(iz,ix,1) = A(iz,ix,2); A(iz,ix,2) = A(iz,ix,3);
    B(iz,ix,1) = B(iz,ix,2); B(iz,ix,2) = B(iz,ix,3);
    A_1(iz,ix,1) = A_1(iz,ix,2); A_1(iz,ix,2) = A_1(iz,ix,3);
    B_1(iz,ix,1) = B_1(iz,ix,2); B_1(iz,ix,2) = B_1(iz,ix,3);
    A_2(iz,ix,1) = A_2(iz,ix,2); A_2(iz,ix,2) = A_2(iz,ix,3);
    B_2(iz,ix,1) = B_2(iz,ix,2); B_2(iz,ix,2) = B_2(iz,ix,3);
    
end  % time loop ends
