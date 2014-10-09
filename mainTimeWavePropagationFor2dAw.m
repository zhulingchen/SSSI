% MAINTIMEWAVEPROPAGATION simulates Kirchoff migration and reverse time
% migration (RTM) with 2-d acoustic wave in time domain with absorbing
% boundary condition (ABC) called Nonsplit Convolutional-PML (CPML)
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Originally written by
% Stuart Kozola for Large Data in MATLAB: A Seismic Data Processing Case Study
% http://www.mathworks.com/matlabcentral/fileexchange/30585-large-data-in-matlab--a-seismic-data-processing-case-study
% The MathWorks Inc., 03/2011
% Modified by
% Lingchen Zhu (zhulingchen@gmail.com), Entao Liu (liuentao@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

close all;
clear;
clc;
%% MicroSeismic Time reversal Example - Layer Cake Model
% in time domain
% modified by Entao Liu
% support arbitrary / random microseismic and receiver positions
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
load('./modelData/velocityModel.mat'); % velocityModel
[nz, nx] = size(velocityModel);

dx = 10;
dz = 10;
x = (1:nx) * dx;
z = (1:nz) * dz;

nBoundary = 20;

% grids and positions of receiver array
recArrType = 'uniform';
idxRecArrLeft = 1;
idxRecArrRight = nx;
nRecs = nx;
if (strcmpi(recArrType, 'uniform'))
    xRecGrid = (idxRecArrLeft:ceil((idxRecArrRight - idxRecArrLeft + 1)/nRecs):idxRecArrRight);
elseif (strcmpi(recArrType, 'random'))
    xRecGrid = (idxRecArrLeft:idxRecArrRight);
    xRecGrid = sort(xRecGrid(randperm(idxRecArrRight - idxRecArrLeft + 1, nRecs)));
else
    error('Receiver array deployment type error!');
end
xRec = xRecGrid * dx;


%% Create shot gathers
% Use the velocity model to simulate a seismic survey.  The wave equations
% is solved using finite differences with a continuous source function
vmin = min(velocityModel(:));
vmax = max(velocityModel(:));

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
nDiffOrder = 2;

% Define frequency parameter for ricker wavelet
f = 20;


%% Generate shot signals
% grids and positions of shot array
nShots = 1;
zShotGrid = 60;
zShot = zShotGrid * dz;
xShotGrid = 50;
xShot = xShotGrid * dx;
delayTimeGrid = 0;

% generate shot source field
sourceTime = zeros([size(V), nt]);
for is = 1:nShots
    % Ricker wavelet
    wave1dTime = ricker(f, nt, dt);
    % % Sinusoid (high frequency components, causing grid dispersion)
    % wave1dTime = sin(2*pi*f*t + 2*pi/nShots * (is-1));
    % wave1dTime(ceil(1/f/dt)+1:end) = 0;
    wave1dTime = [zeros(1, delayTimeGrid(is)), wave1dTime(1:end-delayTimeGrid(is))];
    sourceTime(zShotGrid(is), xShotGrid(is)+nBoundary, :) = reshape(wave1dTime, 1, 1, nt);
end

% plot the velocity model
hFig = figure;
set(hFig, 'PaperPositionMode', 'auto');
figPos = get(hFig, 'Position');

subplot(2,2,1);
imagesc(x, z, velocityModel)
xlabel('Distance (m)'); ylabel('Depth (m)');
title('Velocity Model');
hold on;
hShotPos = plot(xShot, zShot, 'w*');
hold off;
colormap(seismic);

%% Generate shots and save to file and video

filenameVideo = './videos/MSforward.mp4';
if ~exist(filenameVideo, 'file')
    objVideoModelShots = VideoWriter(filenameVideo, 'MPEG-4');
    open(objVideoModelShots);
end

% generate shot record
tic;
[dataTrue, snapshotTrue] = fwdTimeCpmlFor2dAw(V, sourceTime, nDiffOrder, nBoundary, dz, dx, dt);
dataTrue(setdiff(1:nx, xRecGrid)+nBoundary, :) = 0;
timeForward = toc;
fprintf('Generate Forward Timing Record. elapsed time = %fs\n', timeForward);

filenameDataTrue = sprintf('./modelData/dataTrue.mat');
filenameSnapshotTrue = sprintf('./modelData/snapshotTrue.mat');

if ~exist(filenameDataTrue, 'file')
    save(filenameDataTrue, 'dataTrue', '-v7.3');
end

if ~exist(filenameSnapshotTrue, 'file')
    save(filenameSnapshotTrue, 'snapshotTrue', '-v7.3');
end

dataTrue = dataTrue(nBoundary+1:end-nBoundary,:)';

for it = 1:nt
    % plot shot function
    hSub = subplot(2,2,2);
    for is = 1:nShots
        wave1dTime = reshape(sourceTime(zShotGrid(is), xShotGrid(is)+nBoundary, :), 1, nt);
        wave1dTime = mat2gray(wave1dTime);
        wave1dTime = wave1dTime + (is-1);
        plot([1:nt], wave1dTime, it, wave1dTime(it), 'r*'); hold on;
        if (is < nShots)
           plot([1, nt], [is, is], 'k-'); 
        end
    end
    hold off;
    xlim([1, nt]);
    set(hSub, 'YTick', []);
    xlabel('Time'); ylabel('Amplitude');
    title(sprintf('Input source waveform'));
    
    % plot shot record evolution (true)
    ds = zeros(nt, nx);
    ds(1:it, :) = dataTrue(1:it, :);
    subplot(2,2,3)
    imagesc(x, t, ds)
    xlabel('Distance (m)'), ylabel('Time (s)')
    title('Shot Record (True)')
    caxis([-0.1 0.1])
    
    % plot wave propagation (true)
    subplot(2,2,4)
    imagesc(x, z, snapshotTrue(1:end-nBoundary, nBoundary+1:end-nBoundary, it))
    xlabel('Distance (m)'), ylabel('Depth (m)')
    title(sprintf('Wave Propagation (True) t = %.3f', t(it)));
    hold on;
    hShotPos = plot(xShot, zShot, 'w*');
    hold off;
    caxis([-0.14 1])
    
    if exist('objVideoModelShots', 'var')
        writeVideo(objVideoModelShots, im2frame(hardcopy(hFig, '-dzbuffer', '-r0')));
    end
    drawnow;
end %shot loop

if exist('objVideoModelShots', 'var')
    close(objVideoModelShots);
end

%% Time Reversal

filenameVideo = './videos/MSreverse.mp4';
if ~exist(filenameVideo, 'file')
    objVideoModelReverse = VideoWriter(filenameVideo, 'MPEG-4');
    open(objVideoModelReverse);
end

load('./modelData/dataTrue.mat'); % dataTrue

noisyDataTrue = dataTrue;

for ixr = 1:nRecs
    noisyDataTrue(xRecGrid(ixr)+nBoundary,:) = awgn(dataTrue(xRecGrid(ixr)+nBoundary,:), 10, 'measured');
end

tic;
[~, rtmsnapshot] = rvsTimeCpmlFor2dAw(V, noisyDataTrue, nDiffOrder, nBoundary, dz, dx, dt);
timeRT = toc;
fprintf('Generate Reverse Time Record, elapsed time = %fs\n', timeRT);

filenameRTMSnapshot = sprintf('./modelData/rtmsnapshot.mat');

if ~exist(filenameRTMSnapshot, 'file')
    save(filenameRTMSnapshot, 'rtmsnapshot', '-v7.3');
end

noisyDataTrue = noisyDataTrue(nBoundary+1:end-nBoundary,:)';
delete(subplot(2, 2, 2));

for it = nt:-1:30
    
    subplot(2,2,3)
    imagesc(x, t, noisyDataTrue); hold on;
    plot(x, t(it), 'b-', 'LineWidth', 2); hold off;
    xlabel('Distance (m)'), ylabel('Time (s)')
    title('Shot Record (True)')
    caxis([-0.1 0.1])
    
    subplot(2,2,4)
    imagesc(x, z, rtmsnapshot(1:end-nBoundary, nBoundary+1:end-nBoundary, it))
    xlabel('Distance (m)'), ylabel('Depth (m)')
    title(sprintf('Wave Propagation (True) t = %.3f', t(it)));
    
    if exist('objVideoModelReverse', 'var')
        writeVideo(objVideoModelReverse, im2frame(hardcopy(hFig, '-dzbuffer', '-r0')));
    end
    drawnow
end

if exist('objVideoModelReverse', 'var')
    close(objVideoModelReverse);
end