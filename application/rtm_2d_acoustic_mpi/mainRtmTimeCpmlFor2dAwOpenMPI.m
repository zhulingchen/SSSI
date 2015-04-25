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

EPSILON = 1e-6;


%% Set path
run ../setpath;


%% Read in velocity model data and plot it
load([model_data_path, '/velocityModel.mat']); % velocityModel
[nz, nx] = size(velocityModel);

% smooth velocity model using average filter
% filterSmooth = fspecial('average', 5);
% velocityModelSmooth = imfilter(velocityModel, filterSmooth, 'replicate');
load([model_data_path, '/velocityModelSmooth.mat']); % velocityModelSmooth

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
VS = extBoundary(velocityModelSmooth, nBoundary, 2);

% number of approximation order for differentiator operator
nDiffOrder = 5;

% Define frequency parameter for ricker wavelet
f = 20;

%% Set up positions of shot array
idxShotArrLeft = 1;
idxShotArrRight = 100;
nShots = 3;
xShotGrid = (idxShotArrLeft:ceil((idxShotArrRight - idxShotArrLeft + 1)/nShots):idxShotArrRight);

%% Generate shot signals
zShotGrid = 1;
zShot = zShotGrid * dz;
xShot = xShotGrid * dx;

%% Perform RTM in global
hFigRtmProcess = figure;
rtmStacked = zeros(nz+nBoundary, nx+2*nBoundary);
for ixs = 1:nShots
    tic;
    % generate shot source field
    sourceTime = zeros([size(V), nt]);
    
    % Ricker wavelet
    wave1dTime = ricker(f, nt, dt);
    sourceTime(zShotGrid, xShotGrid(ixs)+nBoundary, :) = reshape(wave1dTime, 1, 1, nt);
    
    % Generate the forward snapshot
    [dataTrue, fwdSnapshotTrue] = fwdTimeCpmlFor2dAw(V, sourceTime, nDiffOrder, nBoundary, dz, dx, dt);
    [dataSmooth, fwdSnapshotSmooth] = fwdTimeCpmlFor2dAw(VS, sourceTime, nDiffOrder, nBoundary, dz, dx, dt);
    
    % Acquire difference data for reverse propagation
    dataDelta = dataTrue - dataSmooth;
    
    % Generate the reverse snapshot
    [model, rvsSnapshotDelta] = rvsTimeCpmlFor2dAw(V, dataDelta, nDiffOrder, nBoundary, dz, dx, dt);
    % test begin
    % [model, rvsSnapshotTrue] = rvsTimeCpmlFor2dAw(V, dataTrue, nDiffOrder, nBoundary, dz, dx, dt);
    % test end
    
    % RTM by cross-correlation
    rtmNum = zeros(nz+nBoundary, nx+2*nBoundary);
    rtmDen = EPSILON * ones(nz+nBoundary, nx+2*nBoundary);
    for it = 1:nt
        rtmNum = fwdSnapshotSmooth(:, :, it) .* rvsSnapshotDelta(:, :, it) + rtmNum;
        rtmDen = fwdSnapshotSmooth(:, :, it).^2 + rtmDen;
    end
    rtmStacked = rtmStacked + rtmNum ./ rtmDen;
    
    tRtmPerShot = toc;
    fprintf('Elapsed time of Shot No. %d @ (%d, %d) for RTM = %fs\n', ixs, zShotGrid, xShotGrid(ixs), tRtmPerShot);
    
    figure(hFigRtmProcess);
    imagesc(x, z, rtmStacked(1:end-nBoundary, nBoundary+1:end-nBoundary));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(sprintf('Stacked Image after Shot No. %d', ixs));
    colormap(seismic);
end

filenameRtmStacked = [model_data_path, '/RtmStacked.mat'];
save(filenameRtmStacked, 'rtmStacked', '-v7.3');

%% Perform RTM in MPI
mpi_init;

rtmStacked_local = 0;
for ixs = 1:nShots
    tic;
    
    % generate shot source field
    sourceTime = zeros([size(V), nt]);
    
    % Ricker wavelet
    wave1dTime = ricker(f, nt, dt);
    sourceTime(zShotGrid, xShotGrid(ixs)+nBoundary, :) = reshape(wave1dTime, 1, 1, nt);
    
    % Generate the forward snapshot
    [taskId, dataTrue_mpi, fwdSnapshotTrue_mpi, fwdSnapshotTrue_local] ...
        = fwdTimeCpmlFor2dAw_openmpi_mex(V, sourceTime, nDiffOrder, nBoundary, dz, dx, dt);
    [taskId, dataSmooth_mpi, fwdSnapshotSmooth_mpi, fwdSnapshotSmooth_local] ...
        = fwdTimeCpmlFor2dAw_openmpi_mex(VS, sourceTime, nDiffOrder, nBoundary, dz, dx, dt);
    
    % Acquire difference data for reverse propagation
    if (taskId == 0)
        dataDelta_mpi = dataTrue_mpi - dataSmooth_mpi;
    else
        dataDelta_mpi = zeros(nx+2*nBoundary, nt);
    end
    
    % Generate the reverse snapshot
    [taskId, model_mpi, rvsSnapshotDelta_mpi, rvsSnapshotDelta_local] ...
        = rvsTimeCpmlFor2dAw_openmpi_mex(V, dataDelta_mpi, nDiffOrder, nBoundary, dz, dx, dt);
    
    % RTM by cross-correlation
    rtmNum_local = 0;
    rtmDen_local = EPSILON;
    for it = 1:nt
        rtmNum_local = rtmNum_local + fwdSnapshotSmooth_local(:, :, it) .* rvsSnapshotDelta_local(:, :, it);
        rtmDen_local = rtmDen_local + fwdSnapshotSmooth_local(:, :, it).^2;
    end
    rtmStacked_local = rtmStacked_local + rtmNum_local ./ rtmDen_local;
    
    tRtmPerShot_local = toc;
    fprintf('Elapsed time of Shot No. %d @ (%d, %d) on Processor No. %d for RTM = %fs\n', ixs, zShotGrid, xShotGrid(ixs), taskId, tRtmPerShot_local);
end

filenameRtmStacked_local = [model_data_path, sprintf('/RtmStacked_proc%d.mat', taskId)];
save(filenameRtmStacked_local, 'rtmStacked_local', '-v7.3');

mpi_finalize;

% result comparison
if (taskId == 0)
    delta = dataTrue - dataTrue_mpi;
    fprintf('Forward: Maximum abs difference of data = %d\n', max(abs(delta(:))));
    delta = fwdSnapshotTrue - fwdSnapshotTrue_mpi;
    fprintf('Forward: Maximum abs difference of snapshot = %d\n', max(abs(delta(:))));
    delta = dataSmooth - dataSmooth_mpi;
    fprintf('Forward: Maximum abs difference of data (Smooth Version) = %d\n', max(abs(delta(:))));
    delta = fwdSnapshotSmooth - fwdSnapshotSmooth_mpi;
    fprintf('Forward: Maximum abs difference of snapshot (Smooth Version) = %d\n', max(abs(delta(:))));
    delta = model - model_mpi;
    fprintf('Reverse: Maximum abs difference of model = %d\n', max(abs(delta(:))));
    delta = rvsSnapshotDelta - rvsSnapshotDelta_mpi;
    fprintf('Reverse: Maximum abs difference of snapshot = %d\n', max(abs(delta(:))));
end

%% End Matlab process
quit;
