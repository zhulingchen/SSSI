% MAINWAVETIMECPMLFOR3DAW simulates 3-d acoustic wave propagation in time
% domain with absorbing boundary condition (ABC) called Nonsplit
% Convolutional-PML (CPML)
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Lingchen Zhu (zhulingchen@gmail.com), Entao Liu (liuentao@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo
% Processing Georgia Institute of Technology

close all;
clear;
clc;


%% Seismic Migration Example - Fault Model
% in time domain
% modified by Lingchen Zhu
% support arbitrary / random shot and receiver positions
% support Nonsplit Convolutional-PML (CPML)
EPSILON = 1e-3;


%% Set path
run ../setpath;


%% Read in velocity model data and plot it
nz = 60;
nx = 30;
ny = 40;
dz = 10;
dx = 10;
dy = 10;

% velocity model
% velocityModel = 2500 * ones(nz, nx, ny);
filenameVelocityModel = [model_data_path, '/cake3DP.mat'];
[pathVelocityModel, nameVelocityModel] = fileparts(filenameVelocityModel);
load(filenameVelocityModel);

% smooth velocity model using average filter
% filterSmooth = fspecial('average', 5);
% velocityModelSmooth = imfilter(velocityModel, filterSmooth, 'replicate');
x = (1:nx) * dx;
y = (1:ny) * dy;
z = (1:nz) * dz;
nBoundary = 5;

%% Create shot gathers
% Use the velocity model to simulate a seismic survey.  The wave equations
% is solved using finite differences with a continuous source function
vmin = min(velocityModel(:));
vmax = max(velocityModel(:));

% calculate time step dt from stability crierion for finite difference
% solution of the wave equation.
dt = 0.5*(min([dx, dy, dz])/vmax/sqrt(3));

% determine time samples nt from wave travelime to depth and back to
% surface
nt = round(sqrt((dx*nx)^2 + (dy*ny)^2 + (dz*nz)^2)*2/vmin/dt + 1);
t  = (0:nt-1).*dt;

% add region around model for applying absorbing boundary conditions
V = extBoundary(velocityModel, nBoundary, 3);

% Define frequency parameter for ricker wavelet
f = 20;


%% Generate shots and save to file and video

% shot position
xs = 15 + nBoundary;
ys = 20 + nBoundary;
zs = 30;

% generate shot signal
source = zeros([size(V), nt]);
% source(1, xs, 1) = 1; % impulse input
wave1dTime = ricker(f, nt, dt);
source(zs, xs, ys, :) = reshape(wave1dTime, 1, 1, 1, nt);

% generate shot record
tic;
nDiffOrder = 3;
[dataTrue, snapshotTrue] = fwdTimeCpmlFor3dAw(V, source, nDiffOrder, nBoundary, dz, dx, dy, dt);
timeForward = toc;
fprintf('Generate Forward Timing Record for Shot at x = %dm, y = %dm, z = %dm, time = %fs\n', (xs-nBoundary) * dx, (ys-nBoundary) * dy, zs * dz, timeForward);

%% output
figure;
% plot velocity model
subplot(2, 2, 1);
slice(x, y, z, permute(velocityModel, [3, 2, 1]), ...
    round(linspace(x(2), x(end-1), 5)), ...
    round(linspace(y(2), y(end-1), 5)), ...
    round(linspace(z(2), z(end-1), 10)));
xlabel('X - Easting (m)');
ylabel('Y - Northing (m)');
zlabel('Z - Depth (m)');
title('Velocity Model');
set(gca, 'ZDir', 'reverse');
shading interp;
colormap(seismic);
hold on;
plot3((xs-nBoundary) * dx, (ys-nBoundary) * dy, zs * dz, 'w*');
hold off;

% slice_x = round(linspace(x(2), x(end-1), 3));
% slice_y = round(linspace(y(2), y(end-1), 3));
% slice_z = round(linspace(z(2), z(end-1), 3));

slice_x = median(x);
slice_y = median(y);
slice_z = median(z);

for it = 1:nt
    % plot source function in time domain
    subplot(2, 2, 2);
    plot(t, wave1dTime); hold on;
    plot(t(it), wave1dTime(it), 'r*'); hold off;
    xlim([t(1), t(end)]);
    xlabel('Time (s)'); ylabel('Amplitude');
    title(sprintf('Input source waveform'));
    
    % plot received data traces
    subplot(2, 2, 3);
    dataDisplay = zeros(nx, ny, nt);
    dataDisplay(:, :, 1:it) = dataTrue(nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, 1:it);
    slice(x, y, t, permute(dataDisplay, [2, 1, 3]), ...
        round(linspace(x(2), x(end-1), 5)), ...
        round(linspace(y(2), y(end-1), 5)), ...
        t);
    xlabel('X - Easting (m)');
    ylabel('Y - Northing (m)');
    zlabel('Time (s)');
    title('Shot Record');
    set(gca, 'ZDir', 'reverse');
    shading interp;
    caxis([-0.1, 0.1]);
    
    % plot wave propagation snapshots
    subplot(2, 2, 4);
    slice(x, y, z, permute(snapshotTrue(1:end-nBoundary, nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, it), [3, 2, 1]), ...
        slice_x, slice_y, slice_z);
    xlabel('X - Easting (m)');
    ylabel('Y - Northing (m)');
    zlabel('Z - Depth (m)');
    title(sprintf('Wave Propagation t = %.3fs', t(it)));
    set(gca, 'ZDir', 'reverse');
    shading interp;
    caxis([-0.14, 1]);
    
    drawnow;
end
