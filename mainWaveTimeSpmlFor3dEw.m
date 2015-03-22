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

load('./modelData/homo3DP.mat');
vp = velocityModel;
load('./modelData/homo3DS.mat');
vs = velocityModel;

% dimension check
if (ndims(vp) ~= ndims(vs))
    error('Dimension of P-wave and S-wave velocity models are not the same!');
end
if (any(size(vp) ~= size(vs)))
    error('Dimension of P-wave and S-wave velocity models are not the same!');
end

[nz, nx, ny] = size(velocityModel);

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
zShot = 1 * dz;
xShot = 5 * dx;
yShot = 5 * dy;

% number of approximation order for differentiator operator
nDiffOrder = 3;

hFig = figure;
set(hFig, 'Position', [200, 200, 1000, 500]);
set(hFig, 'PaperPositionMode', 'auto');

% plot velocity model
subplot(2, 3, 1);
slice(x, y, z, permute(vp, [2, 3, 1]), ...
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


%% *************** Check the Stability Condition ***************

if dt > min([dx, dy, dz])/(norm(dCoef(nDiffOrder, 's'), 1) * sqrt(3) * vpmax)
    error('The temporal discretization does not satisfy the Courant-Friedrichs-Lewy sampling criterion to ensure the stability of the FD code!');
end


%% add region around model (vp and vs) for applying absorbing boundary conditions
nBoundary = 20;

VP = extBoundary(vp, nBoundary, 3);
VS = extBoundary(vs, nBoundary, 3);


%% Generate Ricker wavelet as source function
f = 20;               % peak frequency
sourceTime = zeros([size(VP), nt]);
wave1dTime = ricker(f, nt, dt);
sourceTime(zShot/dz, xShot/dx + nBoundary, yShot/dy + nBoundary, :) = reshape(wave1dTime, 1, 1, 1, nt);


%% Generate shots and save to file and video
tic;
[snapshotVzp, snapshotVxp, snapshotVyp, snapshotVzs, snapshotVxs, snapshotVys] = fwdTimeSpmlFor3dEw(VP, VS, sourceTime, nDiffOrder, nBoundary, dz, dx, dy, dt);
timeForward = toc;
fprintf('Generate Forward Timing Record. elapsed time = %fs\n', timeForward);





