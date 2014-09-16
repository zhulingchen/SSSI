% MAINSPARSITYFWIFREQCPMLFOR2DAW simulates the full wave inversion (FWI)
%
% The FWI in frequency domain is used to solve the following problem:
% Given a smooth but obscure velocity model and the received data on
% surface, the true but unknown velocity model is to be approximated by
% estimating the scatter field during iterations.
%
% System background
% ====================================================================================================
%
% m = m_0 + epsilon * delta_m
% True field u: (m(x)(d^2/dt^2) - Laplacian)u(x, t; xs) = -f(x, t; xs)
% Incident field u_0: (m_0(x)(d^2/dt^2) - Laplacian)u_0(x, t; xs) = -f(x, t; xs)
% u = u_0 + u_sc, u_sc is scattered field
%
% Therefore,
% (m_0(x)(d^2/dt^2) - Laplacian)u_sc(y, t; xs) = -epsilon * delta_m(x) * (d^2/dt^2)u(x, t; xs)
% so that
% u_sc = -epsilon * G(\delta_m * (d^2/dt^2)u)
% u = (I + epsilon*A)^(-1)u_0
%   = u_0 - epsilon * A * u_0 + epsilon^2 * A^2 * u_0 - ...
%   = u_0 + epsilon * u_1 + epsilon^2 * u_2 + ...
%
% Therefore, u_1 = -A * u_0, where
% A(f) = G(\delta_m*(d^2/dt^2)f)
% and u_1 satisfies
%
% (m_0(x)(d^2/dt^2) - Laplacian)u_1(y, t; xs) = -delta_m(x) * (d^2/dt^2)u_0(x, t; xs) in time domain
% (-m_0(x)w^2 - Laplacian)U_1(y, jw; xs) = delta_m(x) * w^2 * U_0(x, jw; xs) in frequency domain
% and the solution of U_1 can be written in a linear form as
% U_1 = L * \delta_m
%
% The cost function is:
% J = 1/2 * \sum\w \sum\xs \sum\xr |U_1(xs, xr, jw) - \delta_D(xs, xr, jw)|^2 + lambda * |\delta_m(x)|^2
%
% ====================================================================================================
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
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com), Entao Liu (liuentao@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

close all;
clear;
clc;
%% Full Wave Inversion Example
% in frequency domain

ALPHA = 0.75;
DELTA = 1e-4;
FREQTHRES = 1;
MAXITER = 20;


%% Data source
addpath(genpath('./modelData'));
addpath(genpath('./src'));
if ~isunix
    rmpath(genpath('./src/CurveLab-2.1.3/fdct_usfft_cpp'));
    rmpath(genpath('./src/CurveLab-2.1.3/fdct_wrapping_cpp'));
    rmpath(genpath('./src/CurveLab-2.1.3/fdct3d'));
end


%% Start a pool of Matlab workers
numCores = feature('numcores');
if isempty(gcp('nocreate')) % checking to see if my pool is already open
    myPool = parpool(numCores);
end


%% Read in velocity model data
load('./modelData/velocityModel.mat'); % velocityModel
[nz, nx] = size(velocityModel);
nBoundary = 20;

% smooth velocity model using average filter
load('./modelData/velocityModelSmooth.mat'); % velocityModelSmooth

% a more smooth velocity model for FWI
VS = extBoundary(velocityModelSmooth, nBoundary, 2);
VS = [repmat(VS(1, :), nBoundary, 1); VS];
nAvgSize = [1, 1];
hImageSmooth = fspecial('average', nAvgSize);
VS = imfilter(VS, hImageSmooth);
velocityModelSmooth = VS(nBoundary+1:end-nBoundary, nBoundary+1:end-nBoundary);

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
VS = extBoundary(velocityModelSmooth, nBoundary, 2);

% dimension of frequency-domain solution
nLength = nz * nx;
nLengthWithBoundary = (nz + nBoundary) * (nx + 2*nBoundary);

% number of approximation order for differentiator operator
nDiffOrder = 1;

% Define analog frequency parameter for ricker wavelet
f = 20;
% f = w(550)/(2*pi);


%% Wavelet transform parameters
nlevels_wavelet = 0;            % Decomposition level, all 0 means wavelet
pfilter_wavelet = 'pkva' ;      % Pyramidal filter
dfilter_wavelet = 'pkva' ;      % Directional filter


%% Curvelet transform parameters
is_real = 1;


%% Contourlet transform parameters
nlevels = [2, 3] ;      % Decomposition level
pfilter = 'pkva' ;      % Pyramidal filter
dfilter = 'pkva' ;      % Directional filter


%% Shot data recording at the surface
% generate shot signal
rw1dTime = zeros(1, nt);
for ifreq = 1:length(f)
    rw1dTime = rw1dTime + ricker(f(ifreq), nt, dt);
end
rw1dFreq = fftshift(fft(rw1dTime, nfft), 2);
% find active frequency set with FFT amplitude larger than the threshold
activeW = find(abs(rw1dFreq) > FREQTHRES);
activeW = activeW(activeW ~= nfft / 2 + 1); % skip f = 0Hz

dataTrueFreq = zeros(nRecs, nShots, nfft);
dataDeltaFreq = zeros(nRecs, nShots, nfft);
% receiver positions on extended velocity model
xr = xRecGrid + nBoundary;

% generate shot record and save them in frequency domain
for ixs = 1:nShots %21:nx+20 % shot loop
    
    curXsPos = xShotGrid(ixs) + nBoundary; % shot position on x
    
    % generate shot signal
    sourceTime = zeros([size(V), nt]);
    sourceTime(1, curXsPos, :) = reshape(rw1dTime, 1, 1, nt);
    
    tic;
    [dataTrue, ~] = fwdTimeCpmlFor2dAw(V, sourceTime, nDiffOrder, nBoundary, dz, dx, dt);
    [dataSmooth, ~] = fwdTimeCpmlFor2dAw(VS, sourceTime, nDiffOrder, nBoundary, dz, dx, dt);
    timeForward = toc;
    fprintf('Generate Forward Timing Record for Shot No. %d at x = %dm, elapsed time = %fs\n', curXsPos-nBoundary, x(curXsPos-nBoundary), timeForward);
    
    dataTrue = dataTrue(xr, :);
    dataSmooth = dataSmooth(xr, :);
    
    dataTrueFreq(:, ixs, :) = fftshift(fft(dataTrue, nfft, 2), 2);
    dataDeltaFreq(:, ixs, :) = fftshift(fft(dataTrue - dataSmooth, nfft, 2), 2);
    
end % end shot loop

% save received surface data
filenameDataTrueFreq = './modelData/dataTrueFreq.mat';
filenameDataDeltaFreq = './modelData/dataDeltaFreq0.mat';

if ~exist(filenameDataTrueFreq, 'file')
    save(filenameDataTrueFreq, 'dataTrueFreq', '-v7.3');
end
if ~exist(filenameDataDeltaFreq, 'file')
    save(filenameDataDeltaFreq, 'dataDeltaFreq', '-v7.3');
end

% clear variables and functions from memory
clear('dataTrueFreq');
clear('dataDeltaFreq');
clear('sourceTime');

%% Full wave inversion (FWI)
% (1/v^2)*(d^2)u(z, x, t)/dt^2 + s(z, x, t) = (d^2)u(z, x, t)/dz^2 + (d^2)u(z, x, t)/dx^2
%                                           |
%                                   (Fourier transform), (d^n)f(t)/dt^n -> ((jw)^n)F(jw)
%                                           |
%                                           V
% (w^2)/(v^2)*U(z, x, jw) + (d^2)U(z, x, jw)/dz^2 + (d^2)U(z, x, jw)/dx^2 = S(z, x, jw)
%
% Green's function is the impulse response of the wave equation.


modelOld = zeros(nz, nx);
modelNew = 1./VS(1:end-nBoundary, nBoundary+1:end-nBoundary).^2;

% shot positions on extended velocity model
xs = xShotGrid + nBoundary;

hFig1 = figure(1);
hFig2 = figure(2);
set(hFig2, 'Position', [100, 100, 600, 300]);

iter = 1;
while(norm(modelNew - modelOld, 'fro') / norm(modelOld, 'fro') > DELTA && iter <= MAXITER)
    
    modelOld = modelNew;
    vmOld = 1./sqrt(modelOld);
    vmOld = extBoundary(vmOld, nBoundary, 2);
    load(filenameDataDeltaFreq);
    
    % plot the velocity model
    figure(hFig1);
    imagesc(x, z, vmOld(1:end-nBoundary, nBoundary+1:end-nBoundary));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Previous Velocity Model');
    colormap(seismic); colorbar; caxis manual; caxis([vmin, vmax]);
    
    % Wavelet decomposition
    waveletCoeff = pdfbdec(modelOld, pfilter_wavelet, dfilter_wavelet, nlevels_wavelet);
    [vecWaveletCoeff, sWavelet] = pdfb2vec(waveletCoeff);
    waveletFunc = @(x, mode) pdfb(x, sWavelet, pfilter_wavelet, dfilter_wavelet, nlevels_wavelet, nz, nx, mode);
    
    % Curvelet decomposition
    curveletCoeff = fdct_wrapping(modelOld, is_real);
    [vecCurveletCoeff, sCurvelet] = curvelet2vec(curveletCoeff);
    fdctFunc = @(x, mode) fdct(x, sCurvelet, is_real, nz, nx, mode);
    
    % Contourlet decomposition
    pdfbCoeff = pdfbdec(modelOld, pfilter, dfilter, nlevels);
    [vecPdfbCoeff, sPdfb] = pdfb2vec(pdfbCoeff);
    pdfbFunc = @(x, mode) pdfb(x, sPdfb, pfilter, dfilter, nlevels, nz, nx, mode);
    
    %     %% debug begin
    %     % Wavelet l1-optimization
    %     % lasso: min ||pdfbFunc(x) - b||_2^2 s.t. ||x||_1 < \tau
    %     b = waveletFunc(vecWaveletCoeff, 1);
    %     tau = norm(vecWaveletCoeff, 1);
    %     opts = spgSetParms('verbosity', 1, 'optTol', 1e-12);
    %     vw_spg = spg_lasso(waveletFunc, b, tau, opts);
    %
    %     waveletCoeff_spg = vec2pdfb(vw_spg, sWavelet);
    %     modelOldWavelet_spg = pdfbrec(waveletCoeff_spg, pfilter_wavelet, dfilter_wavelet);
    %     delta = modelOld - modelOldWavelet_spg;
    %     display(max(abs(delta(:))));
    %
    %     % Curvelet l1-optimization
    %     % lasso: min ||fdctFunc(x) - b||_2^2 s.t. ||x||_1 < \tau
    %     b = fdctFunc(vecCurveletCoeff, 1);
    %     tau = norm(vecCurveletCoeff, 1);
    %     opts = spgSetParms('verbosity', 1, 'optTol', 1e-12);
    %     vc_spg = spg_lasso(fdctFunc, b, tau, opts);
    %
    %     curveletCoeff_spg = vec2curvelet(vc_spg, sCurvelet);
    %     modelOldCurvelet_spg = ifdct_wrapping(curveletCoeff_spg, is_real);
    %     delta = modelOld - modelOldCurvelet_spg;
    %     display(max(abs(delta(:))));
    %
    %     % Contourlet l1-optimization
    %     % lasso: min ||pdfbFunc(x) - b||_2^2 s.t. ||x||_1 < \tau
    %     b = pdfbFunc(vecPdfbCoeff, 1);
    %     tau = norm(vecPdfbCoeff, 1);
    %     opts = spgSetParms('verbosity', 1, 'optTol', 1e-12);
    %     vp_spg = spg_lasso(pdfbFunc, b, tau, opts);
    %
    %     pdfbCoeff_spg = vec2pdfb(vp_spg, sPdfb);
    %     modelOldPdfb_spg = pdfbrec(pdfbCoeff_spg, pfilter, dfilter);
    %     delta = modelOld - modelOldPdfb_spg;
    %     display(max(abs(delta(:))));
    %     %% debug end
    
    %     %% debug start
    %     % Wavelet matrix
    %     [v, s] = pdfb2vec(waveletCoeff);
    %     nVec = length(v);
    %     PhiWavelet = zeros(nz * nx, nVec);
    %     for iv = 1:nVec
    %         vtmp = zeros(nVec, 1);
    %         vtmp(iv) = 1;
    %         imtmp = pdfbrec( vec2pdfb(vtmp, s), pfilter_wavelet, dfilter_wavelet ) ;
    %         PhiWavelet(:, iv) = imtmp(:);
    %     end
    %     modelOld2 = PhiWavelet * v;
    %     modelOld2 = reshape(modelOld2, nz, nx);
    %     delta = modelOld2 - modelOld;
    %     display(max(abs(delta(:))));
    %
    %     % Curvelet matrix
    %     [v, s] = curvelet2vec(curveletCoeff);
    %     nVec = length(v);
    %     PhiCurvelet = zeros(nz * nx, nVec);
    %     for iv = 1:nVec
    %         vtmp = zeros(nVec, 1);
    %         vtmp(iv) = 1;
    %         imtmp = ifdct_wrapping( vec2curvelet(vtmp, s), is_real) ;
    %         PhiCurvelet(:, iv) = imtmp(:);
    %     end
    %     modelOld2 = PhiCurvelet * v;
    %     modelOld2 = reshape(modelOld2, nz, nx);
    %     delta = modelOld2 - modelOld;
    %     display(max(abs(delta(:))));
    %
    %     % Contourlet matrix
    %     [v, s] = pdfb2vec(pdfbCoeff);
    %     nVec = length(v);
    %     PhiPdfb = zeros(nz * nx, nVec);
    %     for iv = 1:nVec
    %         vtmp = zeros(nVec, 1);
    %         vtmp(iv) = 1;
    %         imtmp = pdfbrec( vec2pdfb(vtmp, s), pfilter, dfilter ) ;
    %         PhiPdfb(:, iv) = imtmp(:);
    %     end
    %     modelOld2 = PhiPdfb * v;
    %     modelOld2 = reshape(modelOld2, nz, nx);
    %     delta = modelOld2 - modelOld;
    %     display(max(abs(delta(:))));
    %     %% debug end
    
    
    % generate Green's functions
    greenFreqForShotSet = cell(1, length(activeW));
    greenFreqForRecSet = cell(1, length(activeW));
    for idx_w = 1:length(activeW)
        
        iw = activeW(idx_w);
        
        fprintf('Generate Green''s functions for f(%d) = %fHz ... ', iw, w(iw)/(2*pi));
        tic;
        [A, ~] = freqSolverCpml(vmOld, zeros(nz + nBoundary, nx + 2*nBoundary), w(iw), nDiffOrder, nBoundary, dz, dx);
        
        % Green's function for every shot
        sourceFreq = zeros(nLengthWithBoundary, nShots);
        sourceFreq((xs-1)*(nz+nBoundary)+1, :) = eye(nShots, nShots);
        greenFreqForShot = A \ sourceFreq;
        % remove external boundaries
        greenFreqForShot = reshape(greenFreqForShot, nz + nBoundary, nx + 2*nBoundary, nShots);
        greenFreqForShot = greenFreqForShot(1:end-nBoundary, nBoundary+1:end-nBoundary, :);
        greenFreqForShotSet{idx_w} = reshape(greenFreqForShot, nLength, nShots);
        
        % Green's function for every receiver
        sourceFreq = zeros(nLengthWithBoundary, nRecs);
        sourceFreq((xr-1)*(nz+nBoundary)+1, :) = eye(nRecs, nRecs);
        greenFreqForRec = A \ sourceFreq;
        % remove external boundaries
        greenFreqForRec = reshape(greenFreqForRec, nz + nBoundary, nx + 2*nBoundary, nRecs);
        greenFreqForRec = greenFreqForRec(1:end-nBoundary, nBoundary+1:end-nBoundary, :);
        greenFreqForRecSet{idx_w} = reshape(greenFreqForRec, nLength, nRecs);
        
        timePerFreq = toc;
        fprintf('elapsed time = %fs\n', timePerFreq);
        
    end
    
    filenameGreenFreqForShotSet = sprintf('./modelData/greenFreqForShotSet%d.mat', iter);
    save(filenameGreenFreqForShotSet, 'greenFreqForShotSet', '-v7.3');
    
    filenameGreenFreqForRecSet = sprintf('./modelData/greenFreqForRecSet%d.mat', iter);
    save(filenameGreenFreqForRecSet, 'greenFreqForRecSet', '-v7.3');
    
    %% plot updated model optimized in different domains
    figure(hFig2);
    subplot(3, 3, 1);
    imagesc(x, z, velocityModel);
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Actual Velocity Model');
    colormap(seismic); colorbar;
    caxis manual; caxis([vmin, vmax]);
    
    subplot(3, 3, 2);
    imagesc(x, z, vmOld(1:end-nBoundary, nBoundary+1:end-nBoundary));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Smooth Velocity Model');
    colormap(seismic); colorbar;
    caxis manual; caxis([vmin, vmax]);
    
    %% minimization using PQN toolbox in model domain
    func = @(dm) misfitFuncModel(dm, w(activeW), rw1dFreq(activeW), nShots, dataDeltaFreq(:, :, activeW), greenFreqForShotSet, greenFreqForRecSet);
    lowerBound = -inf(nLength,1); %1/vmax^2*ones(nLength,1) - reshape(modelOld, nLength, 1);
    upperBound = inf(nLength,1); %1/vmin^2*ones(nLength,1) - reshape(modelOld, nLength, 1);
    funProj = @(x) boundProject(x, lowerBound, upperBound);
    options.verbose = 3;
    options.optTol = 1e-8;
    options.SPGoptTol = 1e-25;
    options.SPGiters = 100;
    options.adjustStep = 1;
    options.bbInit = 0;
    options.maxIter = 10;
    
    [dm_pqn_model, value_pqn_model] = minConF_PQN_new(func, zeros(nLength, 1), funProj, options);
    
    % updated model
    modelOld = reshape(modelOld, nLength, 1);
    modelNew = modelOld + dm_pqn_model;
    modelOld = reshape(modelOld, nz, nx);
    modelNew = reshape(modelNew, nz, nx);
    modelNew(modelNew < 1/vmax^2) = 1/vmax^2;
    modelNew(modelNew > 1/vmin^2) = 1/vmin^2;
    vmNew = 1./sqrt(modelNew);
    vmNew = extBoundary(vmNew, nBoundary, 2);
    
    subplot(3, 3, 4);
    imagesc(x, z, vmNew(1:end-nBoundary, nBoundary+1:end-nBoundary));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Updated Velocity Model (updated by PQN toolbox in model domain)');
    colormap(seismic); colorbar;
    caxis manual; caxis([vmin, vmax]);
    
    
    %% minimization using PQN toolbox in Wavelet domain
    func = @(dcoeff) misfitFuncSparse(dcoeff, @(x)waveletFunc(x, 1), @(x)waveletFunc(x, 2), w(activeW), rw1dFreq(activeW), nShots, dataDeltaFreq(:, :, activeW), greenFreqForShotSet, greenFreqForRecSet);
    % lowerBound = -inf(length(vecWaveletCoeff), 1); %1/vmax^2*ones(nLength,1) - reshape(modelOld, nLength, 1);
    % upperBound = inf(length(vecWaveletCoeff), 1); %1/vmin^2*ones(nLength,1) - reshape(modelOld, nLength, 1);
    % funProj = @(x) boundProject(x, lowerBound, upperBound);
    tau = norm(vecWaveletCoeff, 1);
    funProj = @(x) sign(x).*projectRandom2C(abs(x), tau);
    options.verbose = 3;
    options.optTol = 1e-8;
    options.SPGoptTol = 1e-25;
    options.SPGiters = 100;
    options.adjustStep = 1;
    options.bbInit = 0;
    options.maxIter = 10;
    
    [dcoeff_pqn_wavelet, value_pqn_wavelet] = minConF_PQN_new(func, zeros(length(vecWaveletCoeff), 1), funProj, options);
    
    dm_pqn_wavelet = real(waveletFunc(dcoeff_pqn_wavelet, 1));
    
    % updated model
    modelOld = reshape(modelOld, nLength, 1);
    modelNew = modelOld + dm_pqn_wavelet;
    modelOld = reshape(modelOld, nz, nx);
    modelNew = reshape(modelNew, nz, nx);
    modelNew(modelNew < 1/vmax^2) = 1/vmax^2;
    modelNew(modelNew > 1/vmin^2) = 1/vmin^2;
    vmNew = 1./sqrt(modelNew);
    vmNew = extBoundary(vmNew, nBoundary, 2);
    
    subplot(3, 3, 5);
    imagesc(x, z, vmNew(1:end-nBoundary, nBoundary+1:end-nBoundary));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Updated Velocity Model (updated by PQN toolbox in Contourlet domain)');
    colormap(seismic); colorbar;
    caxis manual; caxis([vmin, vmax]);
    
    
    %% minimization using PQN toolbox in Curvelet domain
    func = @(dcoeff) misfitFuncSparse(dcoeff, @(x)fdctFunc(x, 1), @(x)fdctFunc(x, 2), w(activeW), rw1dFreq(activeW), nShots, dataDeltaFreq(:, :, activeW), greenFreqForShotSet, greenFreqForRecSet);
    % lowerBound = -inf(length(vecCurveletCoeff), 1); %1/vmax^2*ones(nLength,1) - reshape(modelOld, nLength, 1);
    % upperBound = inf(length(vecCurveletCoeff), 1); %1/vmin^2*ones(nLength,1) - reshape(modelOld, nLength, 1);
    % funProj = @(x) boundProject(x, lowerBound, upperBound);
    tau = norm(vecCurveletCoeff, 1);
    funProj = @(x) sign(x).*projectRandom2C(abs(x), tau);
    options.verbose = 3;
    options.optTol = 1e-8;
    options.SPGoptTol = 1e-25;
    options.SPGiters = 100;
    options.adjustStep = 1;
    options.bbInit = 0;
    options.maxIter = 10;
    
    [dcoeff_pqn_curvelet, value_pqn_curvelet] = minConF_PQN_new(func, zeros(length(vecCurveletCoeff), 1), funProj, options);
    
    dm_pqn_curvelet = real(fdctFunc(dcoeff_pqn_curvelet, 1));
    
    % updated model
    modelOld = reshape(modelOld, nLength, 1);
    modelNew = modelOld + dm_pqn_curvelet;
    modelOld = reshape(modelOld, nz, nx);
    modelNew = reshape(modelNew, nz, nx);
    modelNew(modelNew < 1/vmax^2) = 1/vmax^2;
    modelNew(modelNew > 1/vmin^2) = 1/vmin^2;
    vmNew = 1./sqrt(modelNew);
    vmNew = extBoundary(vmNew, nBoundary, 2);
    
    subplot(3, 3, 6);
    imagesc(x, z, vmNew(1:end-nBoundary, nBoundary+1:end-nBoundary));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Updated Velocity Model (updated by PQN toolbox in Curvelet domain)');
    colormap(seismic); colorbar;
    caxis manual; caxis([vmin, vmax]);
    
    
    %% minimization using PQN toolbox in Contourlet domain
    func = @(dcoeff) misfitFuncSparse(dcoeff, @(x)pdfbFunc(x, 1), @(x)pdfbFunc(x, 2), w(activeW), rw1dFreq(activeW), nShots, dataDeltaFreq(:, :, activeW), greenFreqForShotSet, greenFreqForRecSet);
    % lowerBound = -inf(length(vecPdfbCoeff), 1); %1/vmax^2*ones(nLength,1) - reshape(modelOld, nLength, 1);
    % upperBound = inf(length(vecPdfbCoeff), 1); %1/vmin^2*ones(nLength,1) - reshape(modelOld, nLength, 1);
    % funProj = @(x) boundProject(x, lowerBound, upperBound);
    tau = norm(vecPdfbCoeff, 1);
    funProj = @(x) sign(x).*projectRandom2C(abs(x), tau);
    options.verbose = 3;
    options.optTol = 1e-8;
    options.SPGoptTol = 1e-25;
    options.SPGiters = 100;
    options.adjustStep = 1;
    options.bbInit = 0;
    options.maxIter = 10;
    
    [dcoeff_pqn_pdfb, value_pqn_pdfb] = minConF_PQN_new(func, zeros(length(vecPdfbCoeff), 1), funProj, options);
    
    dm_pqn_pdfb = real(pdfbFunc(dcoeff_pqn_pdfb, 1));
    
    % updated model
    modelOld = reshape(modelOld, nLength, 1);
    modelNew = modelOld + dm_pqn_pdfb;
    modelOld = reshape(modelOld, nz, nx);
    modelNew = reshape(modelNew, nz, nx);
    modelNew(modelNew < 1/vmax^2) = 1/vmax^2;
    modelNew(modelNew > 1/vmin^2) = 1/vmin^2;
    vmNew = 1./sqrt(modelNew);
    vmNew = extBoundary(vmNew, nBoundary, 2);
    
    subplot(3, 3, 7);
    imagesc(x, z, vmNew(1:end-nBoundary, nBoundary+1:end-nBoundary));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Updated Velocity Model (updated by PQN toolbox in Contourlet domain)');
    colormap(seismic); colorbar;
    caxis manual; caxis([vmin, vmax]);
    
    
    %% updated model
    dm = dm_pqn_pdfb;
    
    modelOld = reshape(modelOld, nLength, 1);
    modelNew = modelOld + dm;
    modelOld = reshape(modelOld, nz, nx);
    modelNew = reshape(modelNew, nz, nx);
    modelNew(modelNew < 1/vmax^2) = 1/vmax^2;
    modelNew(modelNew > 1/vmin^2) = 1/vmin^2;
    vmNew = 1./sqrt(modelNew);
    vmNew = extBoundary(vmNew, nBoundary, 2);
    
    figure(3);
    imagesc(x, z, vmNew(1:end-nBoundary, nBoundary+1:end-nBoundary));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Updated Velocity Model');
    colormap(seismic); colorbar; caxis manual; caxis([vmin, vmax]);
    % save current updated velocity model
    filenameVmNew = sprintf('./modelData/vmNew%d.mat', iter);
    save(filenameVmNew, 'vmNew', 'modelNew', '-v7.3');
    
    % clear variables and functions from memory
    clear('greenFreqForShotSet');
    clear('greenFreqForRecSet');
    clear('dataDeltaFreq');
    % load received surface data
    load(filenameDataTrueFreq);
    
    dataDeltaFreq = zeros(nRecs, nShots, nfft);
    
    % update dataDeltaFreq
    for ixs = 1:nShots %21:nx+20 % shot loop
        
        curXsPos = xShotGrid(ixs) + nBoundary; % shot position on x
        
        % generate shot signal
        sourceTime = zeros([size(V), nt]);
        sourceTime(1, curXsPos, :) = reshape(rw1dTime, 1, 1, nt);
        
        % generate shot record
        tic;
        [dataSmooth, ~] = fwdTimeCpmlFor2dAw(vmNew, sourceTime, nDiffOrder, nBoundary, dz, dx, dt);
        timeForward = toc;
        fprintf('Generate Forward Timing Record for Shot No. %d at x = %dm, elapsed time = %fs\n', curXsPos-nBoundary, x(curXsPos-nBoundary), timeForward);
        
        dataSmooth = dataSmooth(xr, :);
        
        dataDeltaFreq(:, ixs, :) = squeeze(dataTrueFreq(:, ixs, :)) - fftshift(fft(dataSmooth, nfft, 2), 2);
        
    end % end shot loop
    filenameDataDeltaFreq = sprintf('./modelData/dataDeltaFreq%d.mat', iter);
    save(filenameDataDeltaFreq, 'dataDeltaFreq', '-v7.3');
    
    % clear variables and functions from memory
    clear('dataTrueFreq');
    clear('dataDeltaFreq');
    clear('sourceTime');
    
    fprintf('Full-wave inversion iteration no. %d, model norm difference = %.6f\n', ...
        iter, norm(modelNew - modelOld, 'fro') / norm(modelOld, 'fro'));
    
    iter = iter + 1;
    
end

%% Terminate the pool of Matlab workers
delete(gcp('nocreate'));

