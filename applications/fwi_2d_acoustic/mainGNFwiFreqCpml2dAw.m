% MAINGNFWIFREQCPML2DAW simulates the full waveform inversion (FWI) with
% 2-d acoustic wave in frequency domain based on the CPML absorbing
% boundary condition and the Gauss-Newton method.
%
% The FWI in frequency domain is used to solve the following problem:
% Given a smooth but obscure velocity model and the received data on
% surface, the true but unknown velocity model is to be approximated by
% estimating the scatter field during iterations.
%
%
% System background
% ====================================================================================================
%
% m = m_0 + delta_m, delta_m is model perturbation
% PDE for true field u:         (m(x)(d^2/dt^2) - Laplacian)u(x, t; xs) = f(x, t; xs)
% PDE for incident field u_0:   (m_0(x)(d^2/dt^2) - Laplacian)u_0(x, t; xs) = f(x, t; xs)
% u = u_0 + u_sc, u_sc is scattered field
%
% Therefore, we have
% (m_0(x)(d^2/dt^2) - Laplacian)u_sc(y, t; xs) = - delta_m(x) * (d^2/dt^2)u(x, t; xs)
%
% In frequency domain
% (- m_0(x)w^2 - Laplacian)U_sc(y, jw; xs) = w^2 * delta_m(x) * U(x, jw; xs)
% whose solution is
% U_sc(y, jw; xs) = w^2 * \sum_x G_0(y, jw; x) * delta_m(x) * U(x, jw; xs)
% where G_0(y, jw; x) is the Green's function of m_0(x) from source x to
% receiver y
%
% U_0(y, jw; xs) = U(y, jw; xs) - U_sc(y, jw; xs)
%                = U(y, jw; xs) - w^2 * \sum_x G_0(y, jw; x) * delta_m(x) * U(x, jw; xs)
%                = (I + A)U(y, jw; xs)
% where operator (matrix) A is composed by (rows of) [- w^2 * G_0(y, jw; x) * delta_m(x)]
%
% Therefore, we have
% U = U_0 + U_sc = (I + A)^{-1} * U_0 ~= U_0 - A * U_0 = U_0 + U_1 by
% discarding quadratic and higher order approximation of the Taylor
% expansion and U_1 = - A * U_0 is the Born approximation of U_sc, i.e.,
% U_1(y, jw; xs) = w^2 * \sum_x G_0(y, jw; x) * delta_m(x) * U_0(x, jw; xs)
% which is the solution of
% (- m_0(x)w^2 - Laplacian)U_1(y, jw; xs) = w^2 * delta_m(x) * U_0(x, jw; xs)
%
% By expanding U_0(x, jw; xs) = G_0(x, jw; xs) * F(x, jw; xs) and the solution
% U_1 can be written in a linear form
% U_1(y, jw; xs) = w^2 * \sum_x F(x, jw; xs) * G_0(x, jw; xs) * G_0(y, jw; x) * delta_m(x), i.e.,
% U_1(y, jw; xs) = L * delta_m(x) where operator (matrix) L is composed of
% (rows of) [w^2 * F(x, jw; xs) * G_0(x, jw; xs) * G_0(y, jw; x)]
%
% The cost function is:
% J = 1/2 * \sum_w \sum_xs \sum_xr |U_1(xr, jw; xs) - U_sc(xr, jw; xs)|^2
%
% ====================================================================================================
%
%
% Purpose
% ====================================================================================================
%
% To find an optimized delta_m(x) such that J is minimized and update the
% velocity model m
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
DELTA = 1e-5;
FREQ_THRES = 1;
NFREQS_PER_BAND = 20;
MAXITER = 20;  % dm is being optimized inside PQN (or L-BFGS) optimization


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
rw1dTime = zeros(nt, 1);
for ifreq = 1:length(f)
    rw1dTime = rw1dTime + ricker(f(ifreq), nt, dt);
end
rw1dFreq = fftshift(fft(rw1dTime, nfft));
% find active frequency set with FFT amplitude larger than the threshold
activeW = find(abs(rw1dFreq) > FREQ_THRES);
activeW = activeW(activeW > nfft / 2 + 1); % choose f > 0Hz
nBands = floor(length(activeW) / NFREQS_PER_BAND);
% set up frequency bands such that FWI is carried out in each band sequentially
activeW = activeW(round(linspace(1, length(activeW), NFREQS_PER_BAND * nBands)));
activeW = reshape(activeW, NFREQS_PER_BAND, nBands);
nFreqs = numel(activeW);

dataTrueFreq = zeros(nRecs, nShots, nFreqs);

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


%% generate shot record in time domain
% dataTrueFreq = zeros(nRecs, nShots, nFreqs);
% parfor idx_shot = 1:nShots % shot loop
%     curXsPos = xShotGrid(idx_shot) + nBoundary; % shot position on x
%     curZsPos = zShotGrid(idx_shot);             % shot position on z
%     
%     % generate shot signal
%     sourceTime = zeros([size(V), nt]);
%     sourceTime(curZsPos, curXsPos, :) = reshape(rw1dTime, 1, 1, nt);
%     
%     tic;
%     [dataTrue, ~] = fwdTimeCpmlFor2dAw(V, sourceTime, nDiffOrder, nBoundary, dz, dx, dt);
%     timeForward = toc;
%     fprintf('Generate Forward Timing Record for Shot No. %d at z = %d, x = %dm, elapsed time = %fs\n', idx_shot, zShot(idx_shot), xShot(idx_shot), timeForward);
%     
%     dataTrue = dataTrue(xr, :);
%     
%     dataTrueFreq_tmp = fftshift(fft(dataTrue, nfft, 2), 2);
%     dataTrueFreq(:, idx_shot, :) = dataTrueFreq_tmp(:, activeW);
% end


%% generate shot record and save them in frequency domain
parfor idx_w = 1:nFreqs
    
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


%% Full wave inversion (FWI)
% (1/v^2)*(d^2)u(z, x, t)/dt^2  = (d^2)u(z, x, t)/dz^2 + (d^2)u(z, x, t)/dx^2 + s(z, x, t)
%                                           |
%                                   (Fourier transform), (d^n)f(t)/dt^n -> ((jw)^n)F(jw)
%                                           |
%                                           V
% (w^2)/(v^2)*U(z, x, jw) + (d^2)U(z, x, jw)/dz^2 + (d^2)U(z, x, jw)/dx^2 = -S(z, x, jw)
%
% Green's function is the impulse response of the wave equation.

modelOld = zeros(nz + nBoundary, nx + 2*nBoundary);
modelNew = MS;

hFigOld = figure(1);
hFigNew = figure(2);


%% FWI main iteration
for iband = 1:nBands
    iter = 1;
    dataTrueFreqCurBand = dataTrueFreq(:, :, (iband-1)*NFREQS_PER_BAND+1:iband*NFREQS_PER_BAND);
    while(norm(modelNew - modelOld, 'fro') / norm(modelOld, 'fro') > DELTA && iter <= MAXITER)
        
        modelOld = modelNew;
        vmOld = sqrt(1./modelOld);
        
        % plot the velocity model
        figure(hFigOld);
        imagesc(x, z, vmOld(1:end-nBoundary, nBoundary+1:end-nBoundary));
        xlabel('Distance (m)'); ylabel('Depth (m)');
        title('Previous Velocity Model');
        colormap(seismic); colorbar; caxis manual; caxis([vmin, vmax]);
        
        %% update dataDeltaFreq based on the new velocity model
        dataDeltaFreqCurBand = zeros(nRecs, nShots, NFREQS_PER_BAND);
        parfor idx_w = 1:NFREQS_PER_BAND
            
            iw = activeW(idx_w, iband);
            
            fprintf('Generate %d frequency responses at f(%d) = %fHz ... ', nShots, iw, w(iw)/(2*pi));
            tic;
            
            % calculate smooth data for all shots in frequency domain for current frequency
            sourceFreq = zeros(nLengthWithBoundary, nShots);
            sourceFreq((xs-1)*(nz+nBoundary)+zs, :) = rw1dFreq(iw) * eye(nShots, nShots);
            [~, snapshotSmoothFreq] = freqCpmlFor2dAw(modelOld, sourceFreq, w(iw), nDiffOrder, nBoundary, dz, dx);
            % get calculated data on the receivers
            dataDeltaFreqCurBand(:, :, idx_w) = dataTrueFreqCurBand(:, :, idx_w) - snapshotSmoothFreq((xr-1)*(nz+nBoundary)+zr, :);
            
            timePerFreq = toc;
            fprintf('elapsed time = %fs\n', timePerFreq);
            
        end
        
        %% generate Green's functions
        greenFreqForShotSet = cell(1, NFREQS_PER_BAND);
        greenFreqForRecSet = cell(1, NFREQS_PER_BAND);
        parfor idx_w = 1:NFREQS_PER_BAND
            
            iw = activeW(idx_w, iband);
            
            fprintf('Generate %d Green''s functions at f(%d) = %fHz ... ', nShots + nRecs, iw, w(iw)/(2*pi));
            tic;
            
            % Green's function for every shot
            sourceFreq = zeros(nLengthWithBoundary, nShots);
            sourceFreq((xs-1)*(nz+nBoundary)+zs, :) = eye(nShots, nShots);
            [~, greenFreqForShotSet{idx_w}] = freqCpmlFor2dAw(modelOld, sourceFreq, w(iw), nDiffOrder, nBoundary, dz, dx);
            
            % Green's function for every receiver
            sourceFreq = zeros(nLengthWithBoundary, nRecs);
            sourceFreq((xr-1)*(nz+nBoundary)+zr, :) = eye(nRecs, nRecs);
            [~, greenFreqForRecSet{idx_w}] = freqCpmlFor2dAw(modelOld, sourceFreq, w(iw), nDiffOrder, nBoundary, dz, dx);
            
            timePerFreq = toc;
            fprintf('elapsed time = %fs\n', timePerFreq);
            
        end
        
        %% minimization using PQN toolbox in model (physical) domain
        func = @(dm) lsBornApproxMisfit(dm, w(activeW(:, iband)), rw1dFreq(activeW(:, iband)), dataDeltaFreqCurBand, ...
            greenFreqForShotSet, greenFreqForRecSet);
        lowerBound = 1e-8 * ones(nLengthWithBoundary, 1) - reshape(modelOld, nLengthWithBoundary, 1); % 1/vmax^2*ones(nLengthWithBoundary, 1) - reshape(modelOld, nLengthWithBoundary, 1);
        upperBound = +inf(nLengthWithBoundary, 1); % 1/vmin^2*ones(nLengthWithBoundary, 1) - reshape(modelOld, nLengthWithBoundary, 1);
        funProj = @(x) boundProject(x, lowerBound, upperBound);
        options.verbose = 3;
        options.optTol = 1e-10;
        options.SPGoptTol = 1e-10;
        options.SPGiters = 5000;
        options.adjustStep = 1;
        options.bbInit = 0;
        options.maxIter = 20;
        
        [dm_pqn_model, misfit_pqn_model] = minConF_PQN_new(func, zeros(nLengthWithBoundary, 1), funProj, options);
        
        %% updated model
        dm = dm_pqn_model;
        misfit = misfit_pqn_model;
        
        modelOld = reshape(modelOld, nLengthWithBoundary, 1);
        modelNew = modelOld + dm;
        modelOld = reshape(modelOld, nz + nBoundary, nx + 2*nBoundary);
        modelNew = reshape(modelNew, nz + nBoundary, nx + 2*nBoundary);
        modelNew(modelNew < 1/vmax^2) = 1/vmax^2;
        modelNew(modelNew > 1/vmin^2) = 1/vmin^2;
        vmNew = sqrt(1./modelNew);
        
        % plot the velocity model
        figure(hFigNew);
        imagesc(x, z, vmNew(1:end-nBoundary, nBoundary+1:end-nBoundary));
        xlabel('Distance (m)'); ylabel('Depth (m)');
        title('Updated Velocity Model');
        colormap(seismic); colorbar; caxis manual; caxis([vmin, vmax]);
        % save current updated velocity model
        filenameVmNew = [pathVelocityModel, sprintf('/vmNew_fband%d_iter%d.mat', iband, iter)];
        save(filenameVmNew, 'vmNew', 'modelNew', 'dm', 'misfit', '-v7.3');
        
        fprintf('Full-wave inversion iteration no. %d, misfit error = %f, model norm difference = %.6f\n', ...
            iter, misfit, norm(modelNew - modelOld, 'fro') / norm(modelOld, 'fro'));
        
        iter = iter + 1;
        
    end
end

%% Terminate the pool of Matlab workers
delete(gcp('nocreate'));

