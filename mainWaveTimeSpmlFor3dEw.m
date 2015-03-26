% MAINWAVETIMESPMLFOR3DEW simulates 3-d elastic wave propagation in time
% domain with absorbing boundary condition (ABC) called Split Perfectly
% Matched Layer (SPML)
%
% 3D Elastic wave propagtaion problem
% Staggered grid
%
% second order for time grid and forth order for spatial grid using staggered grid
% finite diffrence method.
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

close all;
clear;
clc;

addpath(genpath('./modelData'));
addpath(genpath('./src'));


%% Model Parameters setup
% vp = P-wave velocity, vs = S-wave velocity

load('./modelData/cake3DP.mat');
vp = velocityModel;
load('./modelData/cake3DS.mat');
vs = velocityModel;
% vp = 2500 * ones(30, 30, 30);
% vs = vp / sqrt(2);

% dimension check
if (ndims(vp) ~= ndims(vs))
    error('Dimension of P-wave and S-wave velocity models are not the same!');
end
if (any(size(vp) ~= size(vs)))
    error('Dimension of P-wave and S-wave velocity models are not the same!');
end

[nz, nx, ny] = size(vp);

dz = 10;
dx = 10;
dy = 10;
z = (1:nz) * dz;
x = (1:nx) * dx;
y = (1:ny) * dy;

vpmin = min(vp(:));
vpmax = max(vp(:));
vsmin = min(vs(:));
vsmax = max(vs(:));

dt = 0.5*(min([dx, dy, dz])/vpmax/sqrt(3));
nt = round(sqrt((dx*nx)^2 + (dy*ny)^2 + (dz*nz)^2)*2/vsmin/dt + 1);
t  = (0:nt-1)*dt;

% shot position
zShot = 30 * dz;
xShot = 15 * dx;
yShot = 20 * dy;

% number of approximation order for differentiator operator
nDiffOrder = 3;

hFig = figure;
set(hFig, 'Position', [100, 200, 1250, 500]);
set(hFig, 'PaperPositionMode', 'auto');


%% Check the condition of stability
if dt > min([dx, dy, dz])/(norm(dCoef(nDiffOrder, 's'), 1) * sqrt(3) * vpmax)
    error('The temporal discretization does not satisfy the Courant-Friedrichs-Lewy sampling criterion to ensure the stability of the FD code!');
end


%% Add region around model (vp and vs) for applying absorbing boundary conditions
nBoundary = 5;

VP = extBoundary(vp, nBoundary, 3);
VS = extBoundary(vs, nBoundary, 3);


%% Generate Ricker wavelet as source function
f = 20;               % peak frequency
sourceTime = zeros([size(VP), nt]);
wave1dTime = ricker(f, nt, dt);
sourceTime(zShot/dz, xShot/dx + nBoundary, yShot/dy + nBoundary, :) = reshape(wave1dTime, 1, 1, 1, nt);


%% Generate shots
tic;
[snapshotVzp, snapshotVxp, snapshotVyp, snapshotVzs, snapshotVxs, snapshotVys] = fwdTimeSpmlFor3dEw(VP, VS, sourceTime, nDiffOrder, nBoundary, dz, dx, dy, dt);
timeForward = toc;
fprintf('Generate Forward Timing Record. elapsed time = %fs\n', timeForward);


%% Video output
filenameVideo3dEw = './videos/Wave_3D_Elastic.mp4';
if ~exist(filenameVideo3dEw, 'file')
    objVideo3dEw = VideoWriter(filenameVideo3dEw, 'MPEG-4');
    open(objVideo3dEw);
end

% plot velocity model
subplot(2, 4, 1);
slice(x, y, z, permute(vp, [3, 2, 1]), ...
    round(linspace(x(2), x(end-1), 5)), ...
    round(linspace(y(2), y(end-1), 5)), ...
    round(linspace(z(2), z(end-1), 10)));
xlabel('Easting (m)');
ylabel('Northing (m)');
zlabel('Depth (m)');
title('Velocity Model');
set(gca, 'ZDir', 'reverse');
shading interp;
colormap(seismic);
hold on;
plot3(xShot, yShot, zShot, 'w*');
hold off;

% slice_x = round(linspace(x(2), x(end-1), 3));
% slice_y = round(linspace(y(2), y(end-1), 3));
% slice_z = round(linspace(z(2), z(end-1), 3));

slice_x = median(x);
slice_y = median(y);
slice_z = median(z);

for it = 1:nt
    % plot P-wave propagation snapshots
    subplot(2, 4, 2);
    slice(x, y, z, permute(snapshotVzp(1:end-nBoundary, nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, it), [3, 2, 1]), ...
        slice_x, slice_y, slice_z);
    xlabel('X - Easting (m)');
    ylabel('Y - Northing (m)');
    zlabel('Z - Depth (m)');
    title(sprintf('P-wave (z-axis component), t = %.3fs', t(it)));
    set(gca, 'ZDir', 'reverse');
    shading interp;
    caxis([-0.005, 0.05]);
    
    subplot(2, 4, 3);
    slice(x, y, z, permute(snapshotVxp(1:end-nBoundary, nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, it), [3, 2, 1]), ...
        slice_x, slice_y, slice_z);
    xlabel('X - Easting (m)');
    ylabel('Y - Northing (m)');
    zlabel('Z - Depth (m)');
    title(sprintf('P-wave (x-axis component), t = %.3fs', t(it)));
    set(gca, 'ZDir', 'reverse');
    shading interp;
    caxis([-0.005, 0.05]);
    
    subplot(2, 4, 4);
    slice(x, y, z, permute(snapshotVyp(1:end-nBoundary, nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, it), [3, 2, 1]), ...
        slice_x, slice_y, slice_z);
    xlabel('X - Easting (m)');
    ylabel('Y - Northing (m)');
    zlabel('Z - Depth (m)');
    title(sprintf('P-wave (y-axis component), t = %.3fs', t(it)));
    set(gca, 'ZDir', 'reverse');
    shading interp;
    caxis([-0.005, 0.05]);

%     subplot(2, 4, 2);
%     imagesc(squeeze(snapshotVzp(zShot/dz, nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, it)));
%     xlabel('Distance (m)'); ylabel('Depth (m)');
%     title(sprintf('P-wave (z-axis component), t = %.3f', t(it)));
%     caxis([-0.005, 0.05]);
%     
%     subplot(2, 4, 3);
%     imagesc(squeeze(snapshotVxp(zShot/dz, nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, it)));
%     xlabel('Distance (m)'); ylabel('Depth (m)');
%     title(sprintf('P-wave (x-axis component), t = %.3f', t(it)));
%     caxis([-0.005, 0.05]);
%     
%     subplot(2, 4, 4);
%     imagesc(squeeze(snapshotVyp(zShot/dz, nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, it)));
%     xlabel('Distance (m)'); ylabel('Depth (m)');
%     title(sprintf('P-wave (y-axis component), t = %.3f', t(it)));
%     caxis([-0.005, 0.05]);

    % plot source function in time domain
    subplot(2, 4, 5);
    plot([1:nt], wave1dTime); hold on;
    plot(it, wave1dTime(it), 'r*'); hold off;
    xlim([1, nt]);
    xlabel('Time'); ylabel('Amplitude');
    title(sprintf('Input source waveform'));
    
    % plot S-wave propagation snapshots
    subplot(2, 4, 6);
    slice(x, y, z, permute(snapshotVzs(1:end-nBoundary, nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, it), [3, 2, 1]), ...
        slice_x, slice_y, slice_z);
    xlabel('X - Easting (m)');
    ylabel('Y - Northing (m)');
    zlabel('Z - Depth (m)');
    title(sprintf('S-wave (z-axis component), t = %.3fs', t(it)));
    set(gca, 'ZDir', 'reverse');
    shading interp;
    caxis([-0.005, 0.05]);
    
    subplot(2, 4, 7);
    slice(x, y, z, permute(snapshotVxs(1:end-nBoundary, nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, it), [3, 2, 1]), ...
        slice_x, slice_y, slice_z);
    xlabel('X - Easting (m)');
    ylabel('Y - Northing (m)');
    zlabel('Z - Depth (m)');
    title(sprintf('S-wave (x-axis component), t = %.3fs', t(it)));
    set(gca, 'ZDir', 'reverse');
    shading interp;
    caxis([-0.005, 0.05]);
    
    subplot(2, 4, 8);
    slice(x, y, z, permute(snapshotVys(1:end-nBoundary, nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, it), [3, 2, 1]), ...
        slice_x, slice_y, slice_z);
    xlabel('X - Easting (m)');
    ylabel('Y - Northing (m)');
    zlabel('Z - Depth (m)');
    title(sprintf('S-wave (y-axis component), t = %.3fs', t(it)));
    set(gca, 'ZDir', 'reverse');
    shading interp;
    caxis([-0.005, 0.05]);

%     subplot(2, 4, 6);
%     imagesc(squeeze(snapshotVzs(zShot/dz, nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, it)));
%     xlabel('Distance (m)'); ylabel('Depth (m)');
%     title(sprintf('S-wave (z-axis component), t = %.3f', t(it)));
%     caxis([-0.005, 0.05]);
%     
%     subplot(2, 4, 7);
%     imagesc(squeeze(snapshotVxs(zShot/dz, nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, it)));
%     xlabel('Distance (m)'); ylabel('Depth (m)');
%     title(sprintf('S-wave (x-axis component), t = %.3f', t(it)));
%     caxis([-0.005, 0.05]);
%     
%     subplot(2, 4, 8);
%     imagesc(squeeze(snapshotVys(zShot/dz, nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, it)));
%     xlabel('Distance (m)'); ylabel('Depth (m)');
%     title(sprintf('S-wave (y-axis component), t = %.3f', t(it)));
%     caxis([-0.005, 0.05]);
    
    if exist('objVideo3dEw', 'var')
        writeVideo(objVideo3dEw, im2frame(hardcopy(hFig, '-dzbuffer', '-r0')));
    end
    drawnow;
end

% close the VideoWriter object
if exist('objVideo3dEw', 'var')
    close(objVideo3dEw);
end