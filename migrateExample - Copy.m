close all;
clear;
clc;
%% Seismic Migration with a 20 GB Array

%% Data source
% The data used here is from the 2005 EAGE workshop
% dataLocation = pwd; %'C:\DistCompDemos\LargeDataSeismic\benchmark';
addpath('benchmark')
addpath('fileReader')

%% Read in velocity model data
% SegYFileReader open the SEG Y file and reads the header information.
% Creates an object used to access trace data.
tic
V = SegYFileReader('vel_z6.25m_x12.5m_exact.segy',true,false);
toc

%% Read in all traces to memory
tic
V = V(:,:);
toc

imagesc(V), axis equal tight
colormap(seismic)

%% Velocity model parameters
nz = 1911;          % depth samples
dz = 6.25;          % depth spacing (m)
nx = 5395;          % surface samples
dx = 12.5;          % surface spacing (m)

z = (0:nz-1)*dz;
x = (0:nx-1)*dx;

% Modify colormap to show ocean
cmap = colormap(seismic);
cmap(1,:) = [0 0 0];
colormap(cmap)

% update plot
subplot(2,1,1)
imagesc(x,z,V); %axis equal tight
colorbar('NorthOutside')
xlabel('Suface Distance (m)')
ylabel('Depth (m)')

%% Explain Shot Data Source
% Shot records contained in separate 2 GB files.
nShots = 1348;       % total number of shots
dxr = 12.5;          % x-dir spacing between gathers (m)
nxr = 1201;          % number of receivers along x
dt  = 0.006;         % time sample interval (s)
nt  = 2001;          % number of time samples
ds = 50;             % shot increment (m)

xr = (-(nxr-1):1:0)*dxr;
t  = (0:nt-1)*dt;
xshot = (1:nShots)*ds;

%% Load shot records into array
files = {'shots0001_0200.segy','shots0201_0400.segy',...
    'shots0401_0600.segy','shots0601_0800.segy',...
    'shots0801_1000.segy','shots1001_1200.segy',...
    'shots1201_1348.segy'};
shotsInFile = [200:200:1200 1348];
tic
for i = 1:length(files)
    s(i) = SegYFileReader(files{i},true,false);                                                         %#ok<SAGROW>
end
toc

%% Show a shot record on plot
% Plot middle shot record
currentShot = 674;
[fileInd, colInd] = shotRecordLocator(currentShot,s,shotsInFile);
xs = xr + xshot(currentShot);

% load shot record into memory
tic
shot = s{fileInd}(:,colInd);
toc

% plot it
subplot(2,1,2)
imagesc(xs,t,shot), caxis([-.01 .01])
ax = axis; axis([0 x(end) ax(3:4)])
ylabel('Time (s)')
xlabel('Surface Distance (m)')
title(sprintf('Shot Record: %2.0f',currentShot))

% add white window to velocity plot
xrect = [xs(1)  xs(end) xs(end) xs(1) xs(1)];
subplot(2,1,1), hold on
zrect = [z(end) z(end)  z(1)    z(1)  z(end)];
plot(xrect,zrect,'w-')
hold off

%% Loop through shot records and visualize

%createVelocityShotAnimation   % used to create an animation

%% Pull out "tooth" segment for processing
idx = 2200:3400;
x = x(idx);
nx = length(x);
V = V(:,idx);
imagesc(x,z,V), drawnow

%% Find shots in "tooth" Segment
shotsInV = find(xshot >= x(1) & xshot <= x(end));

%% Travel time data: array-like behavior (custom Object)
% Memory mapped file allows me to have an object in MATLAB that behaves
% like an array but only loads the data when needed.
tt = travelTimeMemmap('travelTime.dat',size(V),36);

%% Plot every 25 samples
close all

tic
for i = 1:25:nx
    imagesc(x,z,tt(:,:,i))
    title(['Shot ',num2str(i)]);
    colormap(cmap)
    drawnow;
end
toc

%% Open matlab pool
%workers = 16;
%matlabpool('open','seismic-cluster',workers,...
%           'FileDependencies',{'travelTimeMemmap.m','migrate.m',...
%                              'shot2RecTime'})

%% Migrate 1 shot record
currentShot = 849;
dt = 0.072;
[fileInd, colInd] = shotRecordLocator(currentShot,s,shotsInFile);
shot = s{fileInd}(:,colInd);
M = migrate(tt,diff(shot,2,1),dt,nz,size(shot,2));

imagesc(x,z,M)
axis equal tight
xlabel('Suface Distance (m)')
ylabel('Depth (m)')

%% Migrate all shot records we care about (over tooth region)

% Set up plot
close all
plotProgress(x,z,V,'velocity')
colormap(cmap)
xaxis = axis;
xaxis = xaxis(1:2);

%% Start processing
Stacked = zeros(nz,nx);
ixs =4; 
tstart = tic;
frame = 0;
for currentShot = shotsInV
    tic
    
    %% Plot current shot
    [fileInd, colInd] = shotRecordLocator(currentShot,s,shotsInFile);
    xs = xr + xshot(currentShot);
    
    inV = find(xs >= x(1) & xs <= x(end));
    xs = xs(inV);
    
    shot = s{fileInd}(:, colInd(inV));
    shot = diff(shot,2,1);
    plotProgress(xs,t,shot,'shot','trace',currentShot,...
        'clim',[-.01 .01],'alim',[xaxis 0 t(end)])
    
    %% migrated shot
    shotInd = size(shot,2);
    m = migrate(tt,shot,dt,nz,shotInd,nx);                                                                                 %#ok<*NASGU>
    
    plotProgress(xs,t,m,'migrated','trace',currentShot,...
        'alim',[xaxis 0 t(end)],'clim',[-Inf Inf])
    
    %% stacked image
    Stacked(:,1:shotInd) = Stacked(:,1:shotInd) + m;
    plotProgress(xs,z,diff(Stacked(:,1:shotInd)),'stacked',...
        'trace',currentShot,'alim',[xaxis 0 z(end)*1e-3],...
        'clim',[-Inf Inf],'time',toc(tstart))
        
    drawnow
    
    toc
    
    %frame = frame + 1;
    % movieFrames(frame) = getframe(gcf);
end
toc(tstart)

%% Save movie frames
%save migration56Animation movieFrames

%% Save as avi
% movie2avi(movieFrames,'migrationAnimation')