% MAINFWD3DAW simulates wave propagation with 3-d acoustic wave in time
% domain with absorbing boundary condition (ABC) called Nonsplit
% Convolutional-PML (CPML)
%
% This matlab source file is free for use in academic research. All rights
% reserved.
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


%% Data source
% This example is derived from Gerard Schuster's
% <http://utam.gg.utah.edu/Inter.LAB1/CH2.lab/lab.mig.pre/lab.html MATLAB
% example> and book <http://www.cambridge.org/gb/knowledge/isbn/item2327397/?site_locale=en_GB
% Seismic Interferometry>
%
addpath(genpath('./modelData'));
addpath(genpath('./src'));


%% Read in velocity model data and plot it
nz = 100;
nx = 10;
ny = 10;
dz = 10;
dx = 10;
dy = 10;

% velocity model
% velocityModel = 2500*ones(nz,nx,ny);
load('./modelData/homo3DP.mat');

% smooth velocity model using average filter
% filterSmooth = fspecial('average', 5);
% velocityModelSmooth = imfilter(velocityModel, filterSmooth, 'replicate');
x = (1:nx) * dx;
z = (1:nz) * dz;
y = (1:ny) * dy;
nBoundary = 20;

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
nt = round((sqrt((dx*nx)^2 + (dy*ny)^2 + (dz*nz)^2)*2/vmin/dt + 1) / 4);
t  = (0:nt-1).*dt;

% add region around model for applying absorbing boundary conditions
V = extBoundary(velocityModel, nBoundary, 3);

% Define frequency parameter for ricker wavelet
f = 20;


%% Generate shots and save to file and video

% shot position
xs = 5+nBoundary;
ys = 5+nBoundary;
zs = 1;

% generate shot signal
source = zeros([size(V), nt]);
% source(1, xs, 1) = 1; % impulse input
rw1d = ricker(f, nt, dt);
for i=1:length(rw1d)
    source(zs, xs, ys, i) = rw1d(i);
end

% generate shot record
tic;
nDiffOrder = 2;
[dataTrue, snapshotTrue] = fwdTimeCpmlFor3dAw(V, source, nDiffOrder, nBoundary, dz, dx, dy, dt);
timeForward = toc;
fprintf('Generate Forward Timing Record for Shot at x = %dm, y = %dm, z = %dm, time = %fs\n', (xs-nBoundary) * dx, (ys-nBoundary) * dy, zs * dz, timeForward);

start_t = 1;

for it = start_t:nt
    % plot velocity model
    subplot(2, 2, 1);
    slice(x, y, z, permute(velocityModel, [2, 3, 1]), ...
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
    plot3((xs-nBoundary) * dx, (ys-nBoundary) * dy, zs * dz, 'w*');
    hold off;
    
    % plot source function in time domain
    subplot(2, 2, 2);
    plot([1:nt], rw1d); hold on;
    plot(it, rw1d(it), 'r*'); hold off;
    xlim([1, nt]);
    xlabel('Time'); ylabel('Amplitude');
    colormap(seismic);
    
    % plot received data traces
    subplot(2, 2, 3);
    dataDisplay = zeros(nx, ny, nt);
    dataDisplay(:, :, 1:it) = dataTrue((1:nx)+nBoundary, (1:ny)+nBoundary, 1:it);
    slice(x, y, t, dataDisplay, ...
        round(linspace(x(2), x(end-1), 5)), ...
        round(linspace(y(2), y(end-1), 5)), ...
        t);
    xlabel('Easting (m)');
    ylabel('Northing (m)');
    zlabel('Time (s)');
    title('Shot Record');
    set(gca, 'ZDir', 'reverse');
    shading interp;
    caxis([-0.1 0.1]);
    
    % plot wave propagation snapshots
    subplot(2, 2, 4);
    slice(x, y, z, permute(snapshotTrue(1:end-nBoundary, nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, it), [2, 3, 1]), ...
        round(linspace(x(2), x(end-1), 5)), ...
        round(linspace(y(2), y(end-1), 5)), ...
        z);
    xlabel('Easting (m)');
    ylabel('Northing (m)');
    zlabel('Depth (m)');
    title(sprintf('Wave Propagation t = %.3fs', t(it)));
    set(gca, 'ZDir', 'reverse');
    shading interp;
    caxis([-0.14 1]);
    
    drawnow;
end %shot loop
