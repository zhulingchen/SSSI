% Numerical simulation of elastic wave equation in 2D nonhomogeneous media
% with absorbing bounary conditions (left, right, and bottom, not the surface),
% and FCT (flux-corrected transport) to
% eliminatge the numerical dispersion.
%
% 2D Elastic wave propagtaion problem
% Staggered grid
%
% second order for time grid and forth order for spatial grid using staggered grid
% finite diffrence method.

% Reference:
% Ke-Yang Chen, Finite-Difference simulation of elastic wave with separation in pure P-
% and S-modes, Journal of Computational Methods in Physics, Vol. 2014, Article ID: 108713
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Entao Liu (liuentao@gmail.com), Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

close all;
clear;
clc;


addpath(genpath('./modelData'));
addpath(genpath('./src'));


%% Model Parameters setup
% vp = P-wave velocity, vs = S-wave velocity

load('./modelData/velocityModelP.mat');
vp = velocityModel;
load('./modelData/velocityModelS.mat');
vs = velocityModel;

% dimension check
if (ndims(vp) ~= ndims(vs))
    error('Dimension of P-wave and S-wave velocity models are not the same!');
end
if (any(size(vp) ~= size(vs)))
    error('Dimension of P-wave and S-wave velocity models are not the same!');
end

[nz, nx] = size(velocityModel);

dz = 10;
dx = 10;
x = (1:nx)*dx;
z = (1:nz)*dz;

vpmin = min(vp(:));
vpmax = max(vp(:));
vsmin = min(vs(:));
vsmax = max(vs(:));

dt = 0.75*(dz/vpmax/sqrt(2));
nt = round(sqrt((dx*nx)^2 + (dz*nz)^2)*2/vpmin/dt + 1);
t  = (0:nt-1).*dt;

% shot position
zShot = 50 * dz;
xShot = 20 * dx;

hFig = figure;
set(hFig, 'Position', [200, 200, 1000, 500]);
set(hFig, 'PaperPositionMode', 'auto');

% plot the velocity model
subplot(2,3,1);
imagesc(x, z, vp);
xlabel('Distance (m)'); ylabel('Depth (m)');
title('Velocity Model');
hold on;
hShotPos = plot(xShot, zShot, 'w*');
hold off;
colormap(seismic);


%% *************** Check the Stability Condition **********

if dt*sqrt(vpmax^2/dx^2+vsmax^2/dz^2)>1
    error('Bad Parameters setup. Finite difference grid is not stable.');
end


%% add region around model (vp and vs) for applying absorbing boundary conditions
nBoundary = 20;

VP = extBoundary(vp, nBoundary, 2);
VS = extBoundary(vs, nBoundary, 2);

[nz, nx] = size(VP);          % update nz, nx for the extended model


%% Generate Ricker wavelet as source function
f = 20;               % peak frequency
sourceTime = zeros([size(VP), nt]);
wave1dTime = ricker(f, nt, dt);
sourceTime(zShot/dz, xShot/dx + nBoundary, :) = reshape(wave1dTime, 1, 1, nt);


%% Absorbing boundary condition (ABC): Split PML (SPML)

ixb = 1:nBoundary;          % index of x outside left boundary
ixb2 = nx-nBoundary+ixb;    % index of x outside right boundary
izb  = 1:nz-nBoundary;      % index of z inside down boundary
izb2 = nz-nBoundary+ixb;    % index of z outside down boundary

xDampLeft = dampPml(repmat(fliplr(ixb), nz, 1), VP(:, ixb), nBoundary);
xDampRight = dampPml(repmat(ixb, nz, 1), VP(:, ixb2), nBoundary);
xDamp = [xDampLeft, zeros(nz, nx-2*nBoundary), xDampRight];

zDampUp = dampPml(repmat(fliplr(ixb.'), 1, nx), VP(1:nBoundary, :), nBoundary);
zDampDown = dampPml(repmat(ixb.', 1, nx), VP(izb2, :), nBoundary);
zDamp = [zeros(nz-nBoundary, nx); zDampDown];

% common parameters
VP2=VP.^2;
VS2=VS.^2;
TOX=dt/dx;                   % dt over dx
TOZ=dt/dz;                   % dt over dz


%% Setup variables
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
dataVxp = zeros(nz,nx,nt);
dataVzp = zeros(nz,nx,nt);
dataVxs = zeros(nz,nx,nt);
dataVzs = zeros(nz,nx,nt);

% common indicies
r=3; C = dCoef(r,'s'); % differential coefficients for order 2*r
iz = 1+r:nz-r;      % interior z
ix = 1+r:nx-r;      % interior x


%% ************* Time Iteration *********************
tic;
for it=1:nt
    % In the reference the first index is referred to z direction and
    % second index is referred to x direction. So swap the xDamp and zDamp
    A_1(iz,ix,2) = (1-0.5.*zDamp(iz,ix).*dt)./(1+0.5.*zDamp(iz,ix).*dt).*A_1(iz,ix,1)...
        + TOX./(1+0.5.*zDamp(iz,ix).*dt).*(C(1)*(Vx(iz+1,ix,2)-Vx(iz,ix,2))+C(2)*(Vx(iz+3,ix,2)-Vx(iz-2,ix,2))...
        + C(3)*(Vx(iz+2,ix,2)-Vx(iz-3,ix,2)));
    
    
    A_2(iz,ix,2) = (1-0.5.*xDamp(iz,ix).*dt)./(1+0.5.*xDamp(iz,ix).*dt).*A_2(iz,ix,1)...
        + TOZ./(1+0.5.*xDamp(iz,ix).*dt).*(C(1)*(Vz(iz,ix,2)-Vz(iz,ix-1,2))+C(2)*(Vz(iz,ix+1,2)-Vz(iz,ix-2,2))...
        + C(3)*(Vz(iz,ix+2,2)-Vz(iz,ix-3,2)));
    
    B_1(iz,ix,2) = (1-0.5.*zDamp(iz,ix).*dt)./(1+0.5.*zDamp(iz,ix).*dt).*B_1(iz,ix,1)...
        - TOX./(1+0.5.*zDamp(iz,ix).*dt).*(C(1)*(Vz(iz,ix,2)-Vz(iz-1,ix,2))+C(2)*(Vz(iz+1,ix,2)-Vz(iz-2,ix,2))...
        + C(3)*(Vz(iz+2,ix,2)-Vz(iz-3,ix,2)));
    
    B_2(iz,ix,2) = (1-0.5.*xDamp(iz,ix).*dt)./(1+0.5.*xDamp(iz,ix).*dt).*B_2(iz,ix,1)...
        + TOZ./(1+0.5.*xDamp(iz,ix).*dt).*(C(1)*(Vx(iz,ix+1,2)-Vx(iz,ix,2))+C(2)*(Vx(iz,ix+2,2)-Vx(iz,ix-1,2))...
        + C(2)*(Vx(iz,ix+3,2)-Vx(iz,ix-2,2)));
    
    A(iz,ix,2) = A_1(iz,ix,2) + A_2(iz,ix,2);
    
    B(iz,ix,2) = B_1(iz,ix,2) + B_2(iz,ix,2);
    
    Vxp(iz,ix,3) = (1-0.5.*zDamp(iz,ix).*dt)./(1+0.5.*zDamp(iz,ix).*dt).*Vxp(iz,ix,2)...
        + VP2(iz,ix).*TOX./(1+0.5.*zDamp(iz,ix).*dt).*(C(1).*(A(iz,ix,2)-A(iz-1,ix,2))+C(2).*(A(iz+1,ix,2)-A(iz-2,ix,2))...
        + C(3).*(A(iz+2,ix,2)-A(iz-3,ix,2)))...
        + sourceTime(iz,ix,it);      % soure term
    
    Vzp(iz,ix,3) = (1-0.5.*xDamp(iz,ix).*dt)./(1+0.5.*xDamp(iz,ix).*dt).*Vzp(iz,ix,2)...
        + 0.25.*(VP2(iz+1,ix)+VP2(iz+1,ix+1)+VP2(iz,ix)+VP2(iz,ix+1))...
        .*TOZ./(1+0.5.*xDamp(iz,ix).*dt).*(C(1)*(A(iz,ix+1,2)-A(iz,ix,2))+C(2)*(A(iz,ix+2,2)-A(iz,ix-1,2))...
        + C(3)*(A(iz,ix+3,2)-A(iz,ix-2,2)) );
    
    Vxs(iz,ix,3) = (1-0.5.*xDamp(iz,ix).*dt)./(1+0.5.*xDamp(iz,ix).*dt).*Vxs(iz,ix,2)...
        + VS2(iz,ix).*TOZ./(1+0.5.*xDamp(iz,ix).*dt).*(C(1)*(B(iz,ix,2)-B(iz,ix-1,2))+C(2)*(B(iz,ix+1,2)-B(iz,ix-2,2))...
        + C(3)*(B(iz,ix+2,2)-B(iz,ix-3,2)))...
        + sourceTime(iz,ix,it);      % source term
    
    
    Vzs(iz,ix,3) = (1-0.5.*zDamp(iz,ix).*dt)./(1+0.5.*zDamp(iz,ix).*dt).*Vzs(iz,ix,2)...
        - 0.25.*(VS2(iz,ix)+VS2(iz,ix-1)+VS2(iz-1,ix)+VS2(iz-1,ix-1))...
        .*TOX./(1+0.5.*zDamp(iz,ix).*dt).*(C(1)*(B(iz+1,ix,2)-B(iz,ix,2))+C(2)*(B(iz+2,ix,2)-B(iz-1,ix,2))...
        + C(3)*(B(iz+3,ix,2)-B(iz-2,ix,2)) );
    
    Vx(iz,ix,3) = Vxp(iz,ix,3) + Vxs(iz,ix,3);
    
    Vz(iz,ix,3) = Vzp(iz,ix,3) + Vzs(iz,ix,3);
    
    % ***********  record snapshot for X component and Z component **************************
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
toc;

%************************** Video Making **************************
filenameVideo2dEw = './videos/ElasticWave.mp4';
if ~exist(filenameVideo2dEw, 'file')
    objVideo2dEw = VideoWriter(filenameVideo2dEw, 'MPEG-4');
    open(objVideo2dEw);
end

for it=1:nt
    subplot(2,3,2);
    imagesc(dataVxp(1:end-nBoundary, (nBoundary+1):(100+nBoundary),it));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(sprintf('P-wave (x-axis component), t = %.3f', t(it)));
    caxis([-0.05 0.5]);
    
    subplot(2,3,3);
    imagesc(dataVzp(1:end-nBoundary, (nBoundary+1):(100+nBoundary),it));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(sprintf('P-wave (z-axis component), t = %.3f', t(it)));
    caxis([-0.05 0.5]);
    
    % plot shot function
    subplot(2,3,4);
    plot([1:nt], wave1dTime); hold on;
    plot(it, wave1dTime(it), 'r*'); hold off;
    xlim([1, nt]);
    xlabel('Time'); ylabel('Amplitude');
    title(sprintf('Input source waveform'));
    
    subplot(2,3,5);
    imagesc(dataVxs(1:end-nBoundary, (nBoundary+1):(100+nBoundary),it));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(sprintf('S-wave (x-axis component), t = %.3f', t(it)));
    caxis([-0.05 0.5]);
    
    subplot(2,3,6);
    imagesc(dataVzs(1:end-nBoundary, (nBoundary+1):(100+nBoundary),it));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(sprintf('S-wave (z-axis component), t = %.3f', t(it)));
    caxis([-0.05 0.5]);
    
    if exist('objVideo2dEw', 'var')
        writeVideo(objVideo2dEw, im2frame(hardcopy(hFig, '-dzbuffer', '-r0')));
    end
    drawnow;
end

if exist('objVideo2dEw', 'var')
    close(objVideo2dEw);
end
%************************************************************