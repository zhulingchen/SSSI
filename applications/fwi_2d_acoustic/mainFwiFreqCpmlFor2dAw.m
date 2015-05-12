% MAINFWIFREQCPMLFOR2DAW simulates the full waveform inversion (FWI) with
% 2-d acoustic wave in frequency domain based on the CPML absorbing
% boundary condition
%
% The FWI in frequency domain is used to solve the following problem: Given
% a smooth but obscure velocity model as a starting point and then minimize
% the least-squares misfit function defined by the differences at the
% receiver positions between the recorded seismic data and the modeled
% seismic data for each source-receiver pair of the seismic survey.
%
%
% System background
% ====================================================================================================
%
% The true velocity model m is approached by iteratively running the FWI on
% the currently estimated velocity model m'
%
% PDE for true field u_obs:                    (m(x)(d^2/dt^2) - Laplacian)u_obs(x, t; xs) = f(x, t; xs)
% PDE for calculated (modeled) field u_cal:    (m'(x)(d^2/dt^2) - Laplacian)u_cal(x, t; xs) = f(x, t; xs)
% u_obs(x, t; xs) is the observed (recorded) data at position x caused by
% the shot at position xs
% u_cal(x, t; xs) is the calculated (modeled) data at position x caused by
% the shot at position xs
%
% In frequency domain
% PDE for true field u_obs:                    (- m(x)w^2 - Laplacian)U_obs(x, jw; xs) = F(x, jw; xs)
% PDE for calculated (modeled) field u_cal:    (- m'(x)w^2 - Laplacian)U_cal(x, jw; xs) = F(x, jw; xs)
% whose solution is
% U_cal(x, jw; xs) = G'(x, jw; xs) * F(x, jw; xs)
% or
% U_cal(x, jw; xs) = A(jw; xs) \ (-F(x, jw; xs))
% where G'(y, jw; x) is the Green's function of m' from source xs to
% receiver x and A(jw; xs) is a discretization matrix (Helmholtz operator)
% in frequency domain that maps 2D stencils into its columns
%
%
% The cost function is:
% J = 1/2 * \sum_w \sum_xs \sum_xr |U_cal(xr, jw; xs) - U_obs(xr, jw; xs)|^2
%
% ====================================================================================================
%
%
% Purpose
% ====================================================================================================
%
% To find an estimate m' of m such that J is minimized
%
% The minimum of the misfit function J is sought in the vicinity of the
% starting model m_0', m_{k+1}' = m_k' + a_k * g_k where g_k is the
% gradient of J for m_k'
%
% ====================================================================================================
%
%
% Method
% ====================================================================================================
%
% J is minimized using quasi-Newton method
%
% ====================================================================================================
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

close all;
clear;
clc;


%% Full Wave Inversion Example
% in frequency domain

ALPHA = 0.75;
DELTA = 1e-4;
FREQTHRES = 2;
MAXITER = 1;    % actually no need to do more than 1 iteration outside PQN (or L-BFGS) optimization iterations


%% Set path
run([fileparts(pwd), '/setpath']);


%% Read in velocity model data
filenameVelocityModel = [model_data_path, '/velocityModel.mat'];
[pathVelocityModel, nameVelocityModel] = fileparts(filenameVelocityModel);
load(filenameVelocityModel); % velocityModel
[nz, nx] = size(velocityModel);

% smooth velocity model used average filter
filenameVelocityModelSmooth = [model_data_path, '/velocityModelSmooth.mat'];
load(filenameVelocityModelSmooth); % velocityModelSmooth

nBoundary = 20;

dx = 10;
dz = 10;
x = (1:nx) * dx;
z = (1:nz) * dz;

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
zShotGrid = ones(1, nShots); % shots are on the surface
xShot = xShotGrid * dx;
zShot = zShotGrid * dz;

% grids and positions of receiver array (all on the surface)
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
zRecGrid = ones(1, nRecs); % receivers are on the surface
xRec = xRecGrid * dx;
zRec = zRecGrid * dz;


%% Create shot gathers
% Use the velocity model to simulate a seismic survey.  The wave equations
% is solved using finite differences with a continuous source function
vmin = min(velocityModel(:));
vmax = max(velocityModel(:));

% calculate time step dt from stability crierion for finite difference
% solution of the wave equation.
dt = ALPHA * (dz/vmax/sqrt(2));

% determine time samples nt from wave travelime to depth and back to
% surface
nt = round(sqrt((dx*nx)^2 + (dz*nz)^2)*2/vmin/dt + 1);
t  = (0:nt-1).*dt;
nfft = 2^(nextpow2(nt));
dw = 2*pi/nfft;
w = (-pi:dw:pi-dw)/dt; % analog angular frequency \omega = [-pi, pi)/dt

% add region around model for applying absorbing boundary conditions
V = extBoundary(velocityModel, nBoundary, 2);
M = 1./(V.^2);
VS = extBoundary(velocityModelSmooth, nBoundary, 2);
MS = 1./(VS.^2);

% dimension of frequency-domain solution
nLength = nz * nx;
nLengthWithBoundary = (nz + nBoundary) * (nx + 2*nBoundary);

% number of approximation order for differentiator operator
nDiffOrder = 2;

% Define analog frequency parameter for ricker wavelet
f = 20;


%% Shot data recording at the surface
% generate shot signal
rw1dTime = zeros(1, nt);
for ifreq = 1:length(f)
    rw1dTime = rw1dTime + ricker(f(ifreq), nt, dt);
end
rw1dFreq = fftshift(fft(rw1dTime, nfft), 2);
% find active frequency set with FFT amplitude larger than the threshold
activeW = find(abs(rw1dFreq) > FREQTHRES);
activeW = activeW(activeW > nfft / 2 + 1); % choose f > 0Hz

dataTrueFreq = zeros(nRecs, nShots, length(activeW));

% shot positions on extended velocity model
xs = xShotGrid + nBoundary;
zs = zShotGrid;

% receiver positions on extended velocity model
xr = xRecGrid + nBoundary;
zr = zRecGrid;


%% Start a pool of Matlab workers
numCores = feature('numcores');
if isempty(gcp('nocreate')) % checking to see if my pool is already open
    myPool = parpool(numCores);
end


%% generate shot record and save them in frequency domain
parfor idx_w = 1:length(activeW)
    
    iw = activeW(idx_w);
    
    fprintf('Generate %d frequency responses at f(%d) = %fHz ... ', nShots, iw, w(iw)/(2*pi));
    tic;
    
    % received true data for all shots in frequency domain for current frequency
    sourceFreq = zeros(nLengthWithBoundary, nShots);
    sourceFreq((xs-1)*(nz+nBoundary)+zs, :) = rw1dFreq(iw) * eye(nShots, nShots);
    [~, snapshotTrueFreq] = freqCpmlFor2dAw(M, sourceFreq, w(iw), nDiffOrder, nBoundary, dz, dx);
    % get received data on the receivers
    dataTrueFreq(:, :, idx_w) = snapshotTrueFreq((xr-1)*(nz+nBoundary)+zr, :);
    
    timePerFreq = toc;
    fprintf('elapsed time = %fs\n', timePerFreq);
    
end

% save received surface data
filenameDataTrueFreq = [pathVelocityModel, '/dataTrueFreq.mat'];
save(filenameDataTrueFreq, 'dataTrueFreq', '-v7.3');

% clear variables and functions from memory
clear('dataTrueFreq');


%% Full wave inversion (FWI)
% (1/v^2)*(d^2)u(z, x, t)/dt^2  = (d^2)u(z, x, t)/dz^2 + (d^2)u(z, x, t)/dx^2 + s(z, x, t)
%                                           |
%                                   (Fourier transform), (d^n)f(t)/dt^n -> ((jw)^n)F(jw)
%                                           |
%                                           V
% (w^2)/(v^2)*U(z, x, jw) + (d^2)U(z, x, jw)/dz^2 + (d^2)U(z, x, jw)/dx^2 = -S(z, x, jw)

modelOld = zeros(nz + nBoundary, nx + 2*nBoundary);
modelNew = MS;

hFigOld = figure(1);
hFigNew = figure(2);


%% FWI main iteration
iter = 1;
while(norm(modelNew - modelOld, 'fro') / norm(modelOld, 'fro') > DELTA && iter <= MAXITER)
    
    modelOld = modelNew;
    vmOld = sqrt(1./modelOld);
    load(filenameDataTrueFreq);
    
    % plot the velocity model
    figure(hFigOld);
    imagesc(x, z, vmOld(1:end-nBoundary, nBoundary+1:end-nBoundary));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Previous Velocity Model');
    colormap(seismic); colorbar; caxis manual; caxis([vmin, vmax]);
    
    % test begin
    % [f_opt, g_opt] = lsMisfit(M, w(activeW), rw1dFreq(activeW), dataTrueFreq, nz, nx, xs, zs, xr, zr, nDiffOrder, nBoundary, dz, dx);
    % test end
    
    %% minimization using PQN toolbox in model (physical) domain
    func = @(m) lsMisfit(m, w(activeW), rw1dFreq(activeW), dataTrueFreq, nz, nx, xs, zs, xr, zr, nDiffOrder, nBoundary, dz, dx);
    lowerBound = -inf(nLengthWithBoundary, 1);
    upperBound = +inf(nLengthWithBoundary, 1);
    funProj = @(x) boundProject(x, lowerBound, upperBound);
    options.verbose = 3;
    options.optTol = 1e-10;
    options.SPGoptTol = 1e-10;
    options.SPGiters = 5000;
    options.adjustStep = 1;
    options.bbInit = 0;
    options.maxIter = 100;
    
    [modelNew, misfit_model] = minConF_PQN_new(func, reshape(modelOld, nLengthWithBoundary, 1), funProj, options);
    modelNew = reshape(modelNew, nz + nBoundary, nx + 2*nBoundary);
    vmNew = sqrt(1./modelNew);
    
    % plot the velocity model
    figure(hFigNew);
    imagesc(x, z, vmNew(1:end-nBoundary, nBoundary+1:end-nBoundary));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Updated Velocity Model');
    colormap(seismic); colorbar; caxis manual; caxis([vmin, vmax]);
    % save current updated velocity model
    filenameVmNew = [pathVelocityModel, sprintf('/vmNew%d.mat', iter)];
    save(filenameVmNew, 'vmNew', 'modelNew', '-v7.3');
    
    % clear variables and functions from memory
    clear('dataTrueFreq');
    
    fprintf('Full-wave inversion iteration no. %d, misfit error = %f, model norm difference = %.6f\n', ...
        iter, misfit_model, norm(modelNew - modelOld, 'fro') / norm(modelOld, 'fro'));
    
    iter = iter + 1;
    
end


%% Terminate the pool of Matlab workers
delete(gcp('nocreate'));