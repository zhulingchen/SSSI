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
% test begin
% V = reshape(1:numel(V), size(V));
% sourceTime = reshape(1:prod([size(V), nt]), [size(V), nt]);
% test end
tic; [dataTrue, snapshotTrue] = fwdTimeCpmlFor2dAw(V, sourceTime, nDiffOrder, nBoundary, dz, dx, dt); toc;
% tic; [model, rtmsnapshot] = rvsTimeCpmlFor2dAw(V, dataTrue, nDiffOrder, nBoundary, dz, dx, dt); toc;
tic; [dataTrue_mpi, snapshotTrue_mpi, taskId] = fwdTimeCpmlFor2dAw_openmpi_mex(V, sourceTime, nDiffOrder, nBoundary, dz, dx, dt); toc;
if (taskId == 0)
%     delta = test - test_mpi;
%     display(max(delta, [], 1));
%     fprintf('Maximum abs difference of test = %d\n', max(abs(delta(:))));
    delta = dataTrue - dataTrue_mpi;
    fprintf('Maximum abs difference of data = %d\n', max(abs(delta(:))));
    delta = snapshotTrue - snapshotTrue_mpi;
    fprintf('Maximum abs difference of snapshot = %d\n', max(abs(delta(:))));
%     % test for output
%     save('./test_test.mat', 'test', '-v7.3');
%     save('./test_test_mpi.mat', 'test_mpi', '-v7.3');
%     save('./test_snapshotTrue.mat', 'snapshotTrue', '-v7.3');
%     save('./test_snapshotTrue_mpi.mat', 'snapshotTrue_mpi', '-v7.3');
%     objVideoModelShots = VideoWriter('./test_forward.mp4', 'MPEG-4');
%     open(objVideoModelShots);
%     hFig = figure;
%     for it = 1:nt
%         imagesc(x, z, snapshotTrue_mpi(1:end-nBoundary, nBoundary+1:end-nBoundary, it))
%         caxis([-0.14 1])
%         writeVideo(objVideoModelShots, im2frame(hardcopy(hFig, '-dzbuffer', '-r0')));
%         drawnow;
%     end
%     fprintf('Video output complete!\n');
    % test end
end

%% End Matlab process
quit;
