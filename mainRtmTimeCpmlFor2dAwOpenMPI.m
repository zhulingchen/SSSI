% MAINRTMTIMECPMLFOR2DAWOPENMPI simulates Kirchoff migration and reverse
% time migration (RTM) with 2-d acoustic wave in time domain with absorbing
% boundary condition (ABC) called Nonsplit Convolutional-PML (CPML) using
% OpenMPI to achieve high-performance parallel processing
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

close all;
clear;
clc;

addpath(genpath('./modelData'));
addpath(genpath('./src'));

%% Load velocity model
load('./modelData/velocityModel.mat');
[nz, nx] = size(velocityModel);
vmin = min(velocityModel(:));
vmax = max(velocityModel(:));

dx = 10;
dz = 10;
x = (1:nx) * dx;
z = (1:nz) * dz;
nBoundary = 20;

%% Grids and positions of receiver array
xRecGrid = 1:nx;
xRec = xRecGrid * dx;

%% Simulation parameters
% calculate time step dt from stability crierion for finite difference
% solution of the wave equation.
dt = 0.5*(dz/vmax/sqrt(2));

% determine time samples nt from wave travelime to depth and back to
% surface
nt = round((sqrt((dx*nx)^2 + (dz*nz)^2)*2/vmin/dt + 1));
t  = (0:nt-1).*dt;

% add region around model for applying absorbing boundary conditions
V = extBoundary(velocityModel, nBoundary, 2);

% number of approximation order for differentiator operator
nDiffOrder = 5;

% Define frequency parameter for ricker wavelet
f = 20;

%% Generate shot signals
% shot position
zShotGrid = 21;
zShot = zShotGrid * dz;
xShotGrid = 51;
xShot = xShotGrid * dx;

% generate shot source field
sourceTime = zeros([size(V), nt]);

% Ricker wavelet
wave1dTime = ricker(f, nt, dt);
sourceTime(zShotGrid, xShotGrid+nBoundary, :) = reshape(wave1dTime, 1, 1, nt);

%% Generate the shot record
tic; [dataTrue, snapshotTrue] = fwdTimeCpmlFor2dAw(V, sourceTime, nDiffOrder, nBoundary, dz, dx, dt); toc;
tic; [model, rtmsnapshotTrue] = rvsTimeCpmlFor2dAw(V, dataTrue, nDiffOrder, nBoundary, dz, dx, dt); toc;
% tic; [taskId, dataTrue_mpi, snapshotTrue_mpi, snapshotTrue_local] = fwdTimeCpmlFor2dAw_openmpi_mex(V, sourceTime, nDiffOrder, nBoundary, dz, dx, dt); toc;
tic; [taskId, model_mpi, rtmsnapshotTrue_mpi, rtmsnapshotTrue_local] = rvsTimeCpmlFor2dAw_openmpi_mex(V, dataTrue, nDiffOrder, nBoundary, dz, dx, dt); toc;
% % save each local snapshot video
% objVideoModelShots = VideoWriter(sprintf('./snapshotTrue_local_%d.mp4', taskId), 'MPEG-4');
% open(objVideoModelShots);
% hFig = figure;
% for it = 1:nt
%     imagesc(x, z, snapshotTrue_local(:, :, it));
%     caxis([-0.14 1]);
%     writeVideo(objVideoModelShots, im2frame(hardcopy(hFig, '-dzbuffer', '-r0')));
%     truesize;
%     drawnow;
% end
% fprintf('Video output complete for taskId %d!\n', taskId);
if (taskId == 0)
%     delta = dataTrue - dataTrue_mpi;
%     fprintf('Forward: Maximum abs difference of data = %d\n', max(abs(delta(:))));
%     delta = snapshotTrue - snapshotTrue_mpi;
%     fprintf('Forward: Maximum abs difference of snapshot = %d\n', max(abs(delta(:))));
    delta = model - model_mpi;
    fprintf('Reverse: Maximum abs difference of model = %d\n', max(abs(delta(:))));
    delta = rtmsnapshotTrue - rtmsnapshotTrue_mpi;
    fprintf('Reverse: Maximum abs difference of snapshot = %d\n', max(abs(delta(:))));
end

%% End Matlab process
quit;
