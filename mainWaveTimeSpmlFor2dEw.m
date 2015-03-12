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
x = (1:nx)*dx;
z = (1:nz)*dz;

vpmin = min(vp(:));
vpmax = max(vp(:));
vsmin = min(vs(:));
vsmax = max(vs(:));

dt = 0.75*(dz/vpmax/sqrt(2));
nt = round(sqrt((dx*nx)^2 + (dz*nz)^2)*2/vpmin/dt + 1);
t  = (0:nt-1).*dt;

% shot position
zShot = 1 * dz;
xShot = 1 * dx;

hFig = figure;
set(hFig, 'Position', [200, 200, 1000, 500]);
set(hFig, 'PaperPositionMode', 'auto');

% plot the velocity model
subplot(2,3,1);
imagesc(x, z, vp);
xlabel('Distance (m)'); ylabel('Depth (m)');
title('Velocity Model');
hold on;
hShotPos = plot(xShot, zShot, 'w*');
hold off;
colormap(seismic);


%% *************** Check the Stability Condition **********

if dt*sqrt(vpmax^2/dx^2+vsmax^2/dz^2)>1
    error('Bad Parameters setup. Finite difference grid is not stable.');
end


%% add region around model (vp and vs) for applying absorbing boundary conditions
nBoundary = 20;

VP = extBoundary(vp, nBoundary, 2);
VS = extBoundary(vs, nBoundary, 2);

[nz, nx] = size(VP);          % update nz, nx for the extended model

% number of approximation order for differentiator operator
nDiffOrder = 3;

%% Generate Ricker wavelet as source function
f = 20;               % peak frequency
sourceTime = zeros([size(VP), nt]);
wave1dTime = ricker(f, nt, dt);
sourceTime(zShot/dz, xShot/dx + nBoundary, :) = reshape(wave1dTime, 1, 1, nt);

%% Generate shots and save to file and video
tic;
[dataVzp, dataVxp, dataVzs, dataVxs] = fwdTimeSpmlFor2dEw(VP, VS, sourceTime, nDiffOrder, nBoundary, dz, dx, dt);
timeForward = toc;
fprintf('Generate Forward Timing Record. elapsed time = %fs\n', timeForward);



%% Video output
filenameVideo2dEw = './videos/ElasticWave.mp4';
if ~exist(filenameVideo2dEw, 'file')
    objVideo2dEw = VideoWriter(filenameVideo2dEw, 'MPEG-4');
    open(objVideo2dEw);
end

for it=1:nt
    subplot(2, 3, 2);
    imagesc(dataVxp(1:end-nBoundary, (nBoundary+1):(100+nBoundary),it));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(sprintf('P-wave (x-axis component), t = %.3f', t(it)));
    caxis([-0.05 0.5]);
    
    subplot(2, 3, 3);
    imagesc(dataVzp(1:end-nBoundary, (nBoundary+1):(100+nBoundary),it));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(sprintf('P-wave (z-axis component), t = %.3f', t(it)));
    caxis([-0.05 0.5]);
    
    % plot shot function
    subplot(2, 3, 4);
    plot([1:nt], wave1dTime); hold on;
    plot(it, wave1dTime(it), 'r*'); hold off;
    xlim([1, nt]);
    xlabel('Time'); ylabel('Amplitude');
    title(sprintf('Input source waveform'));
    
    subplot(2, 3, 5);
    imagesc(dataVxs(1:end-nBoundary, (nBoundary+1):(100+nBoundary),it));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(sprintf('S-wave (x-axis component), t = %.3f', t(it)));
    caxis([-0.05 0.5]);
    
    subplot(2, 3, 6);
    imagesc(dataVzs(1:end-nBoundary, (nBoundary+1):(100+nBoundary),it));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(sprintf('S-wave (z-axis component), t = %.3f', t(it)));
    caxis([-0.05 0.5]);
    
    if exist('objVideo2dEw', 'var')
        writeVideo(objVideo2dEw, im2frame(hardcopy(hFig, '-dzbuffer', '-r0')));
    end
    drawnow;
end

if exist('objVideo2dEw', 'var')
    close(objVideo2dEw);
end
%************************************************************