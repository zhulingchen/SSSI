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
load('./modelData/velocityModel.mat'); % velocityModel
[nz, nx] = size(velocityModel);

% smooth velocity model using average filter
% filterSmooth = fspecial('average', 5);
% velocityModelSmooth = imfilter(velocityModel, filterSmooth, 'replicate');
load('./modelData/velocityModelSmooth.mat'); % velocityModelSmooth

dx = 10;
dz = 10;
x = (1:nx) * dx;
z = (1:nz) * dz;

nBoundary = 20;

% grids and positions of shot array
shotArrType = 'uniform';
idxShotArrLeft = 1;
idxShotArrRight = nx;
nShots = nx;
if (strcmpi(shotArrType, 'uniform'))
    xShotGrid = (idxShotArrLeft:ceil((idxShotArrRight - idxShotArrLeft + 1)/nShots):idxShotArrRight);
elseif (strcmpi(shotArrType, 'random'))
    xShotGrid = (idxShotArrLeft:idxShotArrRight);
    xShotGrid = sort(xShotGrid(randperm(idxShotArrRight - idxShotArrLeft + 1, nShots)));
else
    error('Shot array deployment type error!');
end
xShot = xShotGrid * dx;

shotWatchList = [1, ceil(nShots/2), nShots];

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

xShotAndRecGrid = union(xShotGrid, xRecGrid);
nShotsAndRecs = length(xShotAndRecGrid);

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
hShotPos = plot(xShot(1), z(1), 'w*');
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
nt = round(sqrt((dx*nx)^2 + (dz*nz)^2)*2/vmin/dt + 1);
t  = (0:nt-1).*dt;

% add region around model for applying absorbing boundary conditions
V = extBoundary(velocityModel, nBoundary, 2);
VS = extBoundary(velocityModelSmooth, nBoundary, 2);

% number of approximation order for differentiator operator
nDiffOrder = 3;

% Define frequency parameter for ricker wavelet
f = 20;


%% Generate shots and save to file and video

figure(hFig);
colormap(seismic); %bone

filenameVideoModelShots = './videos/ModelShots.mp4';
if ~exist(filenameVideoModelShots, 'file')
    objVideoModelShots = VideoWriter(filenameVideoModelShots, 'MPEG-4');
    open(objVideoModelShots);
end

% generate shot signal
rw1dTime = zeros(1, nt);
for ifreq = 1:length(f)
    rw1dTime = rw1dTime + ricker(f(ifreq), nt, dt);
end

for ixs = 1:nShots %21:nx+20 % shot loop
    
    xs = xShotGrid(ixs) + nBoundary; % shot position on x
    
    % % initial wavefield
    % rw = ricker(f, nx + 2*nBoundary, dt, dt*xs, 0);
    % rw = rw(1:nz+nBoundary, :);
    
    % generate shot signal
    source = zeros([size(V), nt]);
    % source(1, xs, 1) = 1; % impulse input
    source(1, xs, :) = reshape(rw1dTime, 1, 1, nt);
    
    % plot shot position
    set(hShotPos, 'XData', x(xs-nBoundary), 'YData', z(1));
    
    % generate shot record
    tic;
    [dataTrue, snapshotTrue] = fwdTimeCpmlFor2dAw(V, source, nDiffOrder, nBoundary, dz, dx, dt);
    [dataSmooth, snapshotSmooth] = fwdTimeCpmlFor2dAw(VS, source, nDiffOrder, nBoundary, dz, dx, dt);
    timeForward = toc;
    fprintf('Generate Forward Timing Record for Shot No. %d at x = %dm, elapsed time = %fs\n', xs-nBoundary, x(xs-nBoundary), timeForward);
    
    filenameDataTrue = sprintf('./modelData/dataTrue%d.mat', xs-nBoundary);
    filenameDataSmooth = sprintf('./modelData/dataSmooth%d.mat', xs-nBoundary);
    filenameSnapshotSmooth = sprintf('./modelData/snapshotSmooth%d.mat', xs-nBoundary);
    
    if ~exist(filenameDataTrue, 'file')
        save(filenameDataTrue, 'dataTrue', '-v7.3');
    end
    
    if ~exist(filenameDataSmooth, 'file')
        save(filenameDataSmooth, 'dataSmooth', '-v7.3');
    end
    
    if ~exist(filenameSnapshotSmooth, 'file')
        save(filenameSnapshotSmooth, 'snapshotSmooth', '-v7.3');
    end
    
    dataTrue = dataTrue(nBoundary+1:end-nBoundary,:)';
    dataSmooth = dataSmooth(nBoundary+1:end-nBoundary,:)';
    
    if ismember(ixs, shotWatchList)
        start_t = 1;
    else
        start_t = nt;
    end
    
    % start = 1;
    
    for it = start_t:nt
        % plot shot function
        subplot(2,2,2);
        plot([1:nt], rw1dTime); hold on;
        plot(it, rw1dTime(it), 'r*'); hold off;
        xlim([1, nt]);
        xlabel('Time'); ylabel('Amplitude');
        title(sprintf('Shot No. %d at x = %dm', ixs, x(xs-nBoundary)));
        
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
        caxis([-0.14 1])
        
        if exist('objVideoModelShots', 'var')
            writeVideo(objVideoModelShots, im2frame(hardcopy(hFig, '-dzbuffer', '-r0')));
        end
        drawnow;
    end %shot loop
end

if exist('objVideoModelShots', 'var')
    close(objVideoModelShots);
end


%% Traveltime by 2D ray-tracing
% Generate the traveltime field for all z = 0 locations

figure(hFig);
colormap(seismic); %bone

filenameVideoModelTravelTime = './videos/ModelTravelTime.mp4';
if ~exist(filenameVideoModelTravelTime, 'file')
    objVideoModelTravelTime = VideoWriter(filenameVideoModelTravelTime, 'MPEG-4');
    open(objVideoModelTravelTime);
end
travelTime = zeros(nz, nx, nShotsAndRecs);
subplot(2,2,2)
for ixs = 1:nShotsAndRecs
    xs = xShotAndRecGrid(ixs); % shot position on x
    travelTime(:, :, ixs) = ray2d(velocityModel, [1 xs], dx); % [1 xs]: shot position on [z, x]
    imagesc(x, z, travelTime(:, :, ixs));
    xlabel('Distance (m)'), ylabel('Depth (m)');
    title(sprintf('Travel time for Shot No. %d at x = %dm', ixs, x(xs)));
    set(hShotPos, 'XData', x(xs));
    drawnow
    if exist('objVideoModelTravelTime', 'var')
        writeVideo(objVideoModelTravelTime, im2frame(hardcopy(hFig, '-dzbuffer', '-r0')));
    end
end

if exist('objVideoModelTravelTime', 'var')
    close(objVideoModelTravelTime);
end

% save results for later re-use
filenameTravelTime = './modelData/travelTime.mat';
if ~exist(filenameTravelTime, 'file')
    save(filenameTravelTime, 'travelTime', '-v7.3')
end


%% Process Shots - Kirchhoff Migration

figure(hFig);
colormap(seismic); %bone

filenameVideoModelKirchhoff = './videos/ModelKirchhoff.mp4';
if ~exist(filenameVideoModelKirchhoff, 'file')
    objVideoModelKirchhoff = VideoWriter(filenameVideoModelKirchhoff, 'MPEG-4');
    open(objVideoModelKirchhoff);
end
load(filenameTravelTime);
Stacked = zeros(nz, nx);

for ixs = 1:nShots
    xs = xShotGrid(ixs); % shot position on x
    
    load(sprintf('./modelData/dataTrue%d.mat', xs)); % data
    dataTrue = dataTrue(nBoundary+1:end-nBoundary, :)';
    tic;
    M = migrate(travelTime, xRecGrid, xShotAndRecGrid, dataTrue, dt, nz, xs);
    timeKirchoffMigration = toc;
    fprintf('Kirchoff Migration for Shot No. %d at x = %d, elapsed time = %fs\n', xs, x(xs), timeKirchoffMigration);
    
    Stacked = Stacked + M;
    
    subplot(2,2,2)
    imagesc(x, z, Stacked)
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Stacked Image');
    colormap(seismic);
    % caxis([-135 135])
    
    subplot(2,2,3)
    imagesc(x,t,dataTrue)
    xlabel('Distance (m)'); ylabel('Time (s)');
    title(sprintf('Current Shot No. %d at x = %dm', ixs, x(xs)));
    caxis([-0.1 0.1])
    
    subplot(2,2,4)
    imagesc(x,t,M)
    xlabel('Distance (m)'); ylabel('Time (s)');
    title(sprintf('Current Migrated Shot No. %d at x = %dm', ixs, x(xs)));
    colormap(seismic);
    % caxis([-5 5])
    
    set(hShotPos, 'XData', x(xs));
    
    drawnow
    if exist('objVideoModelKirchhoff', 'var')
        writeVideo(objVideoModelKirchhoff, im2frame(hardcopy(hFig, '-dzbuffer', '-r0')));
    end
end

if exist('objVideoModelKirchhoff', 'var')
    close(objVideoModelKirchhoff);
end

filenameKirchoffMigStackedMat = './modelData/KirchoffMigStacked.mat';
if ~exist(filenameKirchoffMigStackedMat, 'file')
    StackedKirchoffMig = diff(Stacked(1:end-nBoundary, nBoundary+1:end-nBoundary), 2, 1);
    save(filenameKirchoffMigStackedMat, 'StackedKirchoffMig', '-v7.3');
end

%% Process Shots - Reverse Time Migration

figure(hFig);
colormap(seismic); %bone

filenameVideoModelRTM = './videos/ModelRTM.mp4';
if ~exist(filenameVideoModelRTM, 'file')
    objVideoModelRTM = VideoWriter(filenameVideoModelRTM, 'MPEG-4');
    open(objVideoModelRTM);
end
Stacked = zeros(nz+nBoundary, nx+2*nBoundary);

for ixs = 1:nShots      %1:nx
    xs = xShotGrid(ixs); % shot position on x
    
    load(sprintf('./modelData/dataTrue%d.mat', xs)); % dataTrue
    load(sprintf('./modelData/dataSmooth%d.mat', xs)); % dataSmooth
    dataDelta = dataTrue - dataSmooth;
    dataTrue = dataTrue(nBoundary+1:end-nBoundary, :)';
    
    tic;
    [~, rtmsnapshot] = rvsTimeCpmlFor2dAw(V, dataDelta, nDiffOrder, nBoundary, dz, dx, dt);
    timeRT = toc;
    fprintf('Generate Reverse Time Record for Shot No. %d at x = %d, elapsed time = %fs\n', xs, x(xs), timeRT);
    
    filenameRTMSnapshot = sprintf('./modelData/rtmsnapshot%d.mat', xs);
    
    if ~exist(filenameRTMSnapshot, 'file')
        save(filenameRTMSnapshot,'rtmsnapshot', '-v7.3');
    end
    
    load(sprintf('./modelData/snapshotSmooth%d.mat', xs)); % snapshotSmooth
    
    M = 0;
    s2 = EPSILON;
    for it = 1:nt
        % M = snapshot(:,:,i).*rtmsnapshot(:,:,nt-i+1) + M;
        M = snapshotSmooth(:, :, it) .* rtmsnapshot(:, :, it) + M;
        s2 = snapshotSmooth(:, :, it).^2 + s2;
        
        if ismember(ixs, shotWatchList)
            
            subplot(2,2,2)
            % imagesc(x, z, diff(M(1:end-nBoundary, nBoundary+1:end-nBoundary) ./ s2(1:end-nBoundary, nBoundary+1:end-nBoundary), 2, 1));
            imagesc(x, z, M(1:end-nBoundary, nBoundary+1:end-nBoundary) ./ s2(1:end-nBoundary, nBoundary+1:end-nBoundary));
            xlabel('Distance (m)'); ylabel('Depth (m)');
            title(sprintf('Current Migrated Shot No. %d at x = %dm', ixs, x(xs)));
            colormap(seismic);
            % caxis([-.05 .05])
            
            %             % plot shot record evolution
            %             subplot(2,2,3)
            %             imagesc(x, t, dataForRTM(nBoundary+1:end-nBoundary, :).'); hold on;
            %             plot([x(1), x(end)], [t(nt-i+1), t(nt-i+1)], 'k', 'LineWidth', 2); hold off;
            %             xlabel('Distance (m)'), ylabel('Time (s)')
            %             title('Shot Record')
            %             caxis([-0.1 0.1])
            
            subplot(2,2,3)
            imagesc(x, z, snapshotSmooth(1:end-nBoundary, nBoundary+1:end-nBoundary, it));
            xlabel('Distance (m)'); ylabel('Depth (m)');
            title(sprintf('Forward Time Wave Propagation (Smooth)\nt = %.3f', t(it)));
            caxis([-0.14 1])
            
            subplot(2,2,4)
            % imagesc(x, z, rtmsnapshot(1:end-nBoundary, nBoundary+1:end-nBoundary, it));
            imagesc(x, z, rtmsnapshot(1:end-nBoundary, nBoundary+1:end-nBoundary, nt-it+1));
            xlabel('Distance (m)'); ylabel('Depth (m)');
            % title(sprintf('Reverse Time Wave Propagation t = %.3f', t(it)));
            title(sprintf('Reverse Time Wave Propagation t = %.3f', t(nt-it+1)));
            caxis([-0.14 1])
            
            drawnow
            if exist('objVideoModelRTM', 'var')
                writeVideo(objVideoModelRTM, im2frame(hardcopy(hFig, '-dzbuffer', '-r0')));
            end
        end
    end
    
    % Stacked = Stacked + M;
    Stacked = Stacked + M ./ s2;
    subplot(2,2,2)
    % imagesc(x, z, diff(Stacked(1:end-nBoundary, nBoundary+1:end-nBoundary), 2, 1));
    imagesc(x, z, Stacked(1:end-nBoundary, nBoundary+1:end-nBoundary));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Stacked Image');
    colormap(seismic);
    % caxis([-30 30])
    
    subplot(2,2,3)
    imagesc(x,t,dataTrue)
    xlabel('Distance (m)'); ylabel('Time (s)');
    title(sprintf('Current Shot No. %d at x = %dm', ixs, x(xs)));
    caxis([-0.1 0.1])
    
    subplot(2,2,4)
    % imagesc(x, t, diff(M(1:end-nBoundary, nBoundary+1:end-nBoundary), 2, 1));
    imagesc(x, t, M(1:end-nBoundary, nBoundary+1:end-nBoundary) ./ s2(1:end-nBoundary, nBoundary+1:end-nBoundary));
    xlabel('Distance (m)'); ylabel('Time (s)');
    title(sprintf('Current Migrated Shot No. %d at x = %dm', ixs, x(xs)));
    colormap(seismic);
    % caxis([-1 1])
    
    set(hShotPos,'XData', x(xs));
    drawnow
    if exist('objVideoModelRTM', 'var')
        writeVideo(objVideoModelRTM, im2frame(hardcopy(hFig, '-dzbuffer', '-r0')));
    end
end

% StackedRTM = diff(Stacked(1:end-nBoundary, nBoundary+1:end-nBoundary), 2, 1);
StackedRTM = Stacked(1:end-nBoundary, nBoundary+1:end-nBoundary);

filenameStackedRtmMat = './modelData/StackedRTM.mat';
if ~exist(filenameStackedRtmMat, 'file')
    save(filenameStackedRtmMat, 'StackedRTM', '-v7.3');
end

if exist('objVideoModelRTM', 'var')
    close(objVideoModelRTM);
end