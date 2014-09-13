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
addpath faultModelData
addpath src


%% Read in velocity model data and plot it
nz = 100;
nx = 100;
ny = 100;
dz = 10;
dx = 10;
dy = 10;

velocityModel= 2500*ones(nz,nx,ny);

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
dt = 0.75*(dz/vmax/sqrt(3));

% determine time samples nt from wave travelime to depth and back to
% surface
nt = round(sqrt((dx*nx)^2 + (dz*nz)^2)*2/vmin/dt + 1);
t  = (0:nt-1).*dt;

% add region around model for applying absorbing boundary conditions
V = extBoundary(velocityModel, nBoundary, 3);

% Define frequency parameter for ricker wavelet
f = 20;


%% Generate shots and save to file and video


for ixs = 1
    
    xs = 50+nBoundary; % shot position
    ys = 50+nBoundary;
    % % initial wavefield
    % rw = ricker(f, nx + 2*nBoundary, dt, dt*xs, 0);
    % rw = rw(1:nz+nBoundary, :);
    
    % generate shot signal
    source = zeros([size(V), nt]);
    % source(1, xs, 1) = 1; % impulse input
    rw1d = ricker(f, nt, dt);
    for i=1:length(rw1d)
        source(1, xs, ys, i) = rw1d(i);
    end
    
    
    % generate shot record
    tic;
    nDiffOrder = 1;
    [dataTrue, snapshotTrue] = fwdTimeCpmlFor3dAw(V, source, nDiffOrder, nBoundary, dz, dx, dy, dt);
    timeForward = toc;
    fprintf('Generate Forward Timing Record for Shot No. %d at x = %dm, time = %f\n', xs-nBoundary, x(xs-nBoundary), timeForward);
    
    
    
    start_t = 1;
    
    for it = start_t:nt
        plot shot function
        subplot(2,2,1);
        imagesc(x, z, snapshotTrue(1:end-nBoundary, nBoundary+1:end-nBoundary, 70, it))
        xlabel('Distance (m)'), ylabel('Depth (m)')
        title(sprintf('Wave Propagation (True) t = %.3f', t(it)));
        caxis([-0.14 1])
        
        subplot(2,2,2);
        imagesc(x, z, snapshotTrue(1:end-nBoundary, 70, nBoundary+1:end-nBoundary, it))
        xlabel('Distance (m)'), ylabel('Depth (m)')
        title(sprintf('Wave Propagation (True) t = %.3f', t(it)));
        caxis([-0.14 1])
        
        
        % plot wave propagation (true)
        subplot(2,2,4)
        imagesc(x, z, snapshotTrue(50, nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary, it))
        xlabel('Distance (m)'), ylabel('Depth (m)')
        title(sprintf('Wave Propagation (True) t = %.3f', t(it)));
        caxis([-0.14 1])
        
        drawnow;
    end %shot loop
end



