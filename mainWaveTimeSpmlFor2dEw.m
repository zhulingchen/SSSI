% MAINWAVETIMESPMLFOR2DEW simulates 2-d elastic wave propagation in time
% domain with absorbing boundary condition (ABC) called Split Perfectly
% Matched Layer (SPML)
%
% 2D Elastic wave propagtaion problem
% Staggered grid
%
% second order for time grid and forth order for spatial grid using staggered grid
% finite diffrence method.

% Reference:
% Ke-Yang Chen, Finite-Difference Simulation of Elastic Wave with
% Separation in Pure P- and S-Modes, Journal of Computational Methods in
% Physics, vol. 2014, Article ID 108713, 14 pages, 2014.
% doi:10.1155/2014/108713
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Entao Liu (liuentao@gmail.com), Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

close all;
clear;
clc;

addpath(genpath('./modelData'));
addpath(genpath('./src'));


%% Model Parameters setup
% vp = P-wave velocity, vs = S-wave velocity

load('./modelData/velocityModelP.mat');
vp = velocityModel;
load('./modelData/velocityModelS.mat');
vs = velocityModel;

% dimension check
if (ndims(vp) ~= ndims(vs))
    error('Dimension of P-wave and S-wave velocity models are not the same!');
end
if (any(size(vp) ~= size(vs)))
    error('Dimension of P-wave and S-wave velocity models are not the same!');
end

[nz, nx] = size(velocityModel);

dz = 10;
dx = 10;
z = (1:nz) * dz;
x = (1:nx) * dx;

vpmin = min(vp(:));
vpmax = max(vp(:));
vsmin = min(vs(:));
vsmax = max(vs(:));

dt = 0.5*(min([dx, dz])/vpmax/sqrt(2));
nt = round(sqrt((dx*nx)^2 + (dz*nz)^2)*2/vsmin/dt + 1);
t  = (0:nt-1).*dt;

% shot position
zShot = 50 * dz;
xShot = 20 * dx;

% number of approximation order for differentiator operator
nDiffOrder = 3;

hFig = figure;
set(hFig, 'Position', [200, 200, 1000, 500]);
set(hFig, 'PaperPositionMode', 'auto');


%% Check the condition of stability
if (dt > min([dx, dz])/(norm(dCoef(nDiffOrder, 's'), 1) * sqrt(2) * vpmax))
    error('The temporal discretization does not satisfy the Courant-Friedrichs-Lewy sampling criterion to ensure the stability of the finite difference method!');
end


%% Add region around model (vp and vs) for applying absorbing boundary conditions
nBoundary = 20;

VP = extBoundary(vp, nBoundary, 2);
VS = extBoundary(vs, nBoundary, 2);


%% Generate Ricker wavelet as source function
f = 20;               % peak frequency
sourceTime = zeros([size(VP), nt]);
wave1dTime = ricker(f, nt, dt);
sourceTime(zShot/dz, xShot/dx + nBoundary, :) = reshape(wave1dTime, 1, 1, nt);


%% Generate shots
tic;
[snapshotVzp, snapshotVxp, snapshotVzs, snapshotVxs] = fwdTimeSpmlFor2dEw(VP, VS, sourceTime, nDiffOrder, nBoundary, dz, dx, dt);
timeForward = toc;
fprintf('Generate Forward Timing Record. elapsed time = %fs\n', timeForward);


%% Video output
filenameVideo2dEw = './videos/Wave_2D_Elastic.mp4';
if ~exist(filenameVideo2dEw, 'file')
    objVideo2dEw = VideoWriter(filenameVideo2dEw, 'MPEG-4');
    open(objVideo2dEw);
end

% plot the velocity model
subplot(2,3,1);
imagesc(x, z, vp);
xlabel('Distance (m)'); ylabel('Depth (m)');
title('Velocity Model');
hold on;
hShotPos = plot(xShot, zShot, 'w*');
hold off;
colormap(seismic);

for it = 1:nt
    % plot P-wave propagation snapshots
    subplot(2, 3, 2);
    imagesc(snapshotVzp(1:end-nBoundary, nBoundary+1:end-nBoundary, it));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(sprintf('P-wave (z-axis component), t = %.3f', t(it)));
    caxis([-2e-5, 1e-4]);
    
    subplot(2, 3, 3);
    imagesc(snapshotVxp(1:end-nBoundary, nBoundary+1:end-nBoundary, it));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(sprintf('P-wave (x-axis component), t = %.3f', t(it)));
    caxis([-2e-5, 1e-4]);
    
    % plot shot function
    subplot(2, 3, 4);
    plot(t, wave1dTime); hold on;
    plot(t(it), wave1dTime(it), 'r*'); hold off;
    xlim([t(1), t(end)]);
    xlabel('Time (s)'); ylabel('Amplitude');
    title(sprintf('Input source waveform'));
    
    % plot S-wave propagation snapshots
    subplot(2, 3, 5);
    imagesc(snapshotVzs(1:end-nBoundary, nBoundary+1:end-nBoundary, it));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(sprintf('S-wave (z-axis component), t = %.3f', t(it)));
    caxis([-2e-5, 1e-4]);
    
    subplot(2, 3, 6);
    imagesc(snapshotVxs(1:end-nBoundary, nBoundary+1:end-nBoundary, it));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(sprintf('S-wave (x-axis component), t = %.3f', t(it)));
    caxis([-2e-5, 1e-4]);
    
    if exist('objVideo2dEw', 'var')
        writeVideo(objVideo2dEw, im2frame(hardcopy(hFig, '-dzbuffer', '-r0')));
    end
    drawnow;
end

% close the VideoWriter object
if exist('objVideo2dEw', 'var')
    close(objVideo2dEw);
end