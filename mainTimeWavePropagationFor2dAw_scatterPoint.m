% MAINRTMTIMECPMLFOR2DAW simulates Kirchoff migration and reverse time
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
velocityModel = 2500 * ones(200, 200);
[nz, nx] = size(velocityModel);
velocityModel(round(nz/2)-5:round(nz/2)+5, round(nx/2)-5:round(nx/2)+5) = 5000; % scatter point

dx = 10;
dz = 10;
x = (1:nx) * dx;
z = (1:nz) * dz;

nBoundary = 20;

% grids and positions of shot array
shotz = [10];
shotx = [1000];
if (length(shotz) ~= length(shotx))
    error('The length of z-axis coordinate list and x-axis coordinate list should be the same!');
end

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


% plot the velocity model
hFig = figure;
set(hFig, 'PaperPositionMode', 'auto');
figPos = get(hFig, 'Position');
set(hFig, 'Position', figPos + [0, 0, 0, ~mod(figPos(end), 2)]);

subplot(2,2,1);
imagesc(x, z, velocityModel)
xlabel('Distance (m)'); ylabel('Depth (m)');
title('Velocity Model');
hold on;
hShotPos = plot(shotx, shotz, 'w*');
hold off;
colormap(seismic);


%% Create shot gathers
% Use the velocity model to simulate a seismic survey.  The wave equations
% is solved using finite differences with a continuous source function
vmin = min(velocityModel(:));
vmax = max(velocityModel(:));

% calculate time step dt from stability crierion for finite difference
% solution of the wave equation.
dt = 0.75*(dz/vmax/sqrt(2));

% determine time samples nt from wave travelime to depth and back to
% surface
nt = round(0.5*(sqrt((dx*nx)^2 + (dz*nz)^2)*2/vmin/dt + 1));
t  = (0:nt-1).*dt;

% add region around model for applying absorbing boundary conditions
V = extBoundary(velocityModel, nBoundary, 2);


% number of approximation order for differentiator operator
nDiffOrder = 2;

% Define frequency parameter for ricker wavelet
f = 25;


%% Generate shots and save to file and video

figure(hFig);
colormap(seismic); %bone

filenameVideo = './videos/MSforward_scatterPoint.mp4';
if ~exist(filenameVideo, 'file')
    objVideoModelShots = VideoWriter(filenameVideo, 'MPEG-4');
    open(objVideoModelShots);
end

% generate shot signal
rw1dTime = zeros(1, nt);
for ifreq = 1:length(f)
    rw1dTime = rw1dTime + 0.1*ricker(f(ifreq), nt, dt);
end

% generate shot signal
source = zeros([size(V), nt]);
% source(1, xs, 1) = 1; % impulse input
for is = 1:length(shotz)
    source(shotz(is)/dz, shotx(is)/dx+nBoundary, :) = reshape(rw1dTime, 1, 1, nt);
end

% generate shot record
tic;
[dataTrue, snapshotTrue] = fwdTimeCpmlFor2dAw(V, source, nDiffOrder, nBoundary, dz, dx, dt);
timeForward = toc;
fprintf('Generate Forward Timing Record. elapsed time = %fs\n', timeForward);

filenameDataTrue = sprintf('./modelData/scatterPointModelData/dataTrue.mat');
filenameSnapshotTrue = sprintf('./modelData/scatterPointModelData/snapshotTrue.mat');

if ~exist(filenameDataTrue, 'file')
    save(filenameDataTrue, 'dataTrue', '-v7.3');
end

if ~exist(filenameSnapshotTrue, 'file')
    save(filenameSnapshotTrue, 'snapshotTrue', '-v7.3');
end

dataTrue = dataTrue(nBoundary+1:end-nBoundary,:)';

for it = 1:nt
    % plot shot function
    subplot(2,2,2);
    plot([1:nt], rw1dTime); hold on;
    plot(it, rw1dTime(it), 'r*'); hold off;
    xlim([1, nt]);
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
    hShotPos = plot(shotx, shotz, 'w*');
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

filenameVideo = './videos/MSreverse_scatterPoint.mp4';
if ~exist(filenameVideo, 'file')
    objVideoModelReverse = VideoWriter(filenameVideo, 'MPEG-4');
    open(objVideoModelReverse);
end

load('./modelData/scatterPointModelData/dataTrue.mat'); % dataTrue

noisyDataTrue = dataTrue;

for i=1:nRecs
    noisyDataTrue(i+nBoundary,:) = awgn(dataTrue(i+nBoundary,:), 1000, 'measured');
end

tic;
[~, rtmsnapshot] = rvsTimeCpmlFor2dAw(V, noisyDataTrue, nDiffOrder, nBoundary, dz, dx, dt);
timeRT = toc;
fprintf('Generate Reverse Time Record, elapsed time = %fs\n', timeRT);

filenameRTMSnapshot = sprintf('./modelData/scatterPointModelData/rtmsnapshot.mat');

if ~exist(filenameRTMSnapshot, 'file')
    save(filenameRTMSnapshot, 'rtmsnapshot', '-v7.3');
end

noisyDataTrue = noisyDataTrue(nBoundary+1:end-nBoundary,:)';

for it = nt:-1:20

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
    hold on;
    hShotPos = plot(shotx, shotz, 'w*');
    hold off;
    % caxis([-0.14 1])

    if exist('objVideoModelReverse', 'var')
        writeVideo(objVideoModelReverse, im2frame(hardcopy(hFig, '-dzbuffer', '-r0')));
    end
    drawnow
end

if exist('objVideoModelReverse', 'var')
    close(objVideoModelReverse);
end

%% Reverse Time Migration (RTM)

load('./modelData/scatterPointModelData/snapshotTrue.mat');
load('./modelData/scatterPointModelData/rtmsnapshot.mat');

filenameVideo = './videos/MSrtm_scatterPoint.mp4';
if ~exist(filenameVideo, 'file')
    objVideoRtm = VideoWriter(filenameVideo, 'MPEG-4');
    open(objVideoRtm);
end

hFig = figure;
set(hFig, 'PaperPositionMode', 'auto');

M = 0;
s2 = EPSILON;
for it = 1:nt
    M = snapshotTrue(:, :, it) .* rtmsnapshot(:, :, it) + M;
    s2 = snapshotTrue(:, :, it).^2 + s2;
    
    subplot(2, 2, 1);
    imagesc(x, z, snapshotTrue(1:end-nBoundary, nBoundary+1:end-nBoundary, it));
    xlabel('Distance (m)'), ylabel('Depth (m)');
    title('Forward wave propagation');
    colorbar;
    
    subplot(2, 2, 2);
    imagesc(x, z, rtmsnapshot(1:end-nBoundary, nBoundary+1:end-nBoundary, it));
    xlabel('Distance (m)'), ylabel('Depth (m)');
    title('Reverse wave propagation');
    colorbar;
    
    subplot(2, 2, 3);
    imagesc(x, z, snapshotTrue(1:end-nBoundary, nBoundary+1:end-nBoundary, it) .* rtmsnapshot(1:end-nBoundary, nBoundary+1:end-nBoundary, it));
    xlabel('Distance (m)'), ylabel('Depth (m)');
    title('Instant correlation');
    colorbar;
    
    subplot(2, 2, 4);
    imagesc(x, z, M(1:end-nBoundary, nBoundary+1:end-nBoundary) ./ s2(1:end-nBoundary, nBoundary+1:end-nBoundary));
    xlabel('Distance (m)'), ylabel('Depth (m)');
    title('Reverse time migration');
    colorbar;
    
    if exist('objVideoRtm', 'var')
        writeVideo(objVideoRtm, im2frame(hardcopy(hFig, '-dzbuffer', '-r0')));
    end
    drawnow
    
%     if (mod(it, 50) == 0)
%         fprintf('it = %d\n', it);
%     end
end

if exist('objVideoRtm', 'var')
    close(objVideoRtm);
end