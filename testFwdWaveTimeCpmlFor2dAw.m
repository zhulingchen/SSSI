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
nDiffOrder = 3;

% Define frequency parameter for ricker wavelet
f = 20;

%% Generate shot signals
% shot position
zShotGrid = 1;
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
tic; [dataTrue_mpi, snapshotTrue_mpi, taskId] = fwdTimeCpmlFor2dAw_openmpi_mex(V, sourceTime, nDiffOrder, nBoundary, dz, dx, dt); toc;
display(taskId);
display(size(dataTrue_mpi));
display(size(snapshotTrue_mpi));

%% End Matlab process
quit;
