% MAINSEISMICDATADENOISING performs seismic data denoising using sparse
% dictionary learning
%
% Referred Paper:
% Rubinstein, R.; Zibulevsky, M.; Elad, M., "Double Sparsity: Learning
% Sparse Dictionaries for Sparse Signal Approximation," Signal Processing,
% IEEE Transactions on , vol.58, no.3, pp.1553,1564, March 2010
% doi: 10.1109/TSP.2009.2036477
% keywords: {image coding;image denoising;sparse matrices;3D image
% denoising;computed tomography;double sparsity;learning sparse
% dictionaries;signal representation;sparse coding;sparse signal
% approximation;Computed tomography;K-SVD;dictionary learning;signal
% denoising;sparse coding;sparse representation},
% URL:
% http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5325694&isnumber=5410625
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


%% Data source
addpath(genpath('./modelData'));
addpath(genpath('./src'));
if ~isunix
    rmpath(genpath('./src/CurveLab-2.1.3/fdct_usfft_cpp'));
    rmpath(genpath('./src/CurveLab-2.1.3/fdct_wrapping_cpp'));
    rmpath(genpath('./src/CurveLab-2.1.3/fdct3d'));
end

% dataFile = './modelData/denoising/barbara.mat'; % 512 * 512
% dataFile = './modelData/denoising/timodel_shot_data_II_shot001-320.mat'; % 1332 * 656
% dataFile = './modelData/denoising/frst_ch40hz.mat'; % 400 * 400
% dataFile = './modelData/denoising/bp_eage_shot675.mat'; % 2001 * 1201
% dataFile = './modelData/denoising/bp_eage_shot701.mat'; % 2001 * 1201
% dataFile = './modelData/denoising/bppublic_seg_aa_SEG_free_shots.mat'; % 626 * 176
dataFile = './modelData/denoising/bppublic_2_5d_shots.mat'; % 384 * 252
% dataFile = './modelData/denoising/Overthrust_Sam_int_shots.mat'; % 513 * 301
% dataFile = './modelData/denoising/viking_seismic.mat'; % 1500 * 120
[dataFileDir, dataFileName] = fileparts(dataFile);
load(dataFile); % shot data from Hess VTI synthetic datasets

% % try barbara
% load('./modelData/denoising/barbara.mat'); % dataTrue
% dataTrue = dataTrue(1:256, end-255:end);

nBoundary = 20;
% dataTrue = dataTrue(nBoundary+1:end-nBoundary,:)';
dataTrue = dataTrue(1:16*floor(size(dataTrue, 1)/16), 1:16*floor(size(dataTrue, 2)/16));
[nSamples, nRecs] = size(dataTrue);


%% Remove mean and normalization
minData = min(dataTrue, [], 1);
maxData = max(dataTrue, [], 1);
meanData = mean(dataTrue, 1);
% dataTrue = bsxfun(@times, bsxfun(@minus, dataTrue, minData), 1./abs(maxData - minData));
dataTrue = bsxfun(@times, bsxfun(@minus, dataTrue, meanData), 1./abs(maxData - minData));
% dataTrue = dataTrue / max(dataTrue(:));

% %% Normalize data to unit norm
% for ir = 1:nRecs
%     dataTrue(:, ir) = dataTrue(:, ir) / norm(dataTrue(:, ir), 2);
% end


%% Prepare noisy data
% sampling rate
ts = 9.9e-3;
fs = 1/ts;
% generate colored noise
forder = 10;
fp1 = 0;     % first passband frequency
fp2 = 30;    % second passband frequency
% d = fdesign.bandpass('N,Fp1,Fp2,Ap', forder, fp1, fp2, .5, fs);
% Hd = design(d, 'cheby1');
d = fdesign.lowpass('N,Fp,Fst', forder, (2*fp2/1.05)/fs, 2*fp2/fs);
Hd = design(d, 'firls');
sigma = 0.1;
noise = sigma * randn(size(dataTrue));
noise_filtered = filter(Hd, noise);
noisyData = dataTrue + noise_filtered;
save(fullfile(dataFileDir, [dataFileName, '_noisyData.mat']), 'noisyData', '-v7.3');
trainData = noisyData;


%% Plot figures
hFigDataTrue = figure; imagesc(1:nRecs, 1:nSamples, dataTrue); colormap(gray); colorbar; truesize;
xlabel('Trace Index'); ylabel('Time Sample Index');
title('Original Seismic Data');

hFigNoisyData = figure; imagesc(1:nRecs, 1:nSamples, noisyData); colormap(gray); colorbar; truesize;
psnrNoisyData = 20*log10(sqrt(numel(noisyData)) / norm(dataTrue(:) - noisyData(:), 2));
ssimNoisyData = ssim(noisyData, dataTrue);
xlabel('Trace Index'); ylabel('Time Sample Index');
title(sprintf('Noisy Seismic Data, PSNR = %.2fdB', psnrNoisyData));
saveas(hFigNoisyData, fullfile(dataFileDir, [dataFileName, '_noisyData']), 'fig');
fprintf('------------------------------------------------------------\n');
fprintf('Noisy Seismic Data, PSNR = %.2fdB, SSIM = %.4f\n', psnrNoisyData, ssimNoisyData);


%% Reference: denoising using wavelet
nlevels_wavelet = [0, 0, 0];        % Decomposition level, all 0 means wavelet
pfilter_wavelet = '9/7';            % Pyramidal filter
dfilter_wavelet = 'pkva';           % Directional filter
[vecWaveletCoeff, str] = pdfb2vec(pdfbdec(noisyData, pfilter_wavelet, dfilter_wavelet, nlevels_wavelet));

% Thresholding
% noiseVar = pdfb_nest(size(dataTrue, 1), size(dataTrue, 2), pfilter_wavelet, dfilter_wavelet, nlevels_wavelet);
thres_wavelet = 3 * sigma;
vecWaveletCoeff = vecWaveletCoeff .* (abs(vecWaveletCoeff) > thres_wavelet);

% Reconstruction
cleanData_wavelet = pdfbrec(vec2pdfb(vecWaveletCoeff, str), pfilter_wavelet, dfilter_wavelet);
save(fullfile(dataFileDir, [dataFileName, '_cleanData_wavelet.mat']), 'cleanData_wavelet', '-v7.3');
diffData_wavelet = dataTrue - cleanData_wavelet;
save(fullfile(dataFileDir, [dataFileName, '_diffData_wavelet.mat']), 'diffData_wavelet', '-v7.3');

% Plot figures and PSNR output
hFigCleanedDataWavelet = figure; imagesc(1:nRecs, 1:nSamples, cleanData_wavelet); colormap(gray); colorbar; truesize;
psnrCleanData_wavelet = 20*log10(sqrt(numel(cleanData_wavelet)) / norm(dataTrue(:) - cleanData_wavelet(:), 2));
ssimCleanData_wavelet = ssim(cleanData_wavelet, dataTrue);
xlabel('Trace Index'); ylabel('Time Sample Index');
title(sprintf('Denoised Seismic Data (Wavelet), PSNR = %.2fdB', psnrCleanData_wavelet));
saveas(hFigCleanedDataWavelet, fullfile(dataFileDir, [dataFileName, '_cleanData_wavelet']), 'fig');

hFigDiffDataWavelet = figure; imagesc(1:nRecs, 1:nSamples, diffData_wavelet); colormap(gray); colorbar; truesize;
xlabel('Trace Index'); ylabel('Time Sample Index');
saveas(hFigDiffDataWavelet, fullfile(dataFileDir, [dataFileName, '_diffData_wavelet']), 'fig');

fprintf('------------------------------------------------------------\n');
fprintf('Denoised Seismic Data (Wavelet), PSNR = %.2fdB, SSIM = %.4f\n', psnrCleanData_wavelet, ssimCleanData_wavelet);


%% Reference: denoising using Contourlet
nlevels_contourlet = [2, 3, 4];     % Decomposition level, all 0 means wavelet
pfilter_contourlet = '9/7';         % Pyramidal filter
dfilter_contourlet = 'pkva';        % Directional filter
coeffContourlet = pdfbdec(noisyData, pfilter_contourlet, dfilter_contourlet, nlevels_contourlet);
[vecContourletCoeff, str] = pdfb2vec(coeffContourlet);

% Set up thresholds for coarse scales
noiseVar = pdfb_nest(size(dataTrue, 1), size(dataTrue, 2), pfilter_contourlet, dfilter_contourlet, nlevels_contourlet);
thres_contourlet = 3 * sigma * sqrt(noiseVar.');

% Slightly different thresholds for the finest scale
finestScale = str(end, 1);
finestScaleSize = sum(prod(str(str(:, 1) == finestScale, 3:4), 2));
thres_contourlet(end-finestScaleSize+1:end) = (4/3) * thres_contourlet(end-finestScaleSize+1:end);

% Thresholding
vecContourletCoeff = vecContourletCoeff .* (abs(vecContourletCoeff) > thres_contourlet);

% Reconstruction
cleanData_contourlet = pdfbrec(vec2pdfb(vecContourletCoeff, str), pfilter_contourlet, dfilter_contourlet);
save(fullfile(dataFileDir, [dataFileName, '_cleanData_contourlet.mat']), 'cleanData_contourlet', '-v7.3');
diffData_contourlet = dataTrue - cleanData_contourlet;
save(fullfile(dataFileDir, [dataFileName, '_diffData_contourlet.mat']), 'diffData_contourlet', '-v7.3');

% Plot figures and PSNR output
hFigCleanedDataContourlet = figure; imagesc(1:nRecs, 1:nSamples, cleanData_contourlet); colormap(gray); colorbar; truesize;
psnrCleanData_contourlet = 20*log10(sqrt(numel(cleanData_contourlet)) / norm(dataTrue(:) - cleanData_contourlet(:), 2));
ssimCleanData_contourlet = ssim(cleanData_contourlet, dataTrue);
xlabel('Trace Index'); ylabel('Time Sample Index');
title(sprintf('Denoised Seismic Data (Contourlet, Thresholding), PSNR = %.2fdB', psnrCleanData_contourlet));
saveas(hFigCleanedDataContourlet, fullfile(dataFileDir, [dataFileName, '_cleanData_contourlet']), 'fig');

hFigDiffDataContourlet = figure; imagesc(1:nRecs, 1:nSamples, diffData_contourlet); colormap(gray); colorbar; truesize;
xlabel('Trace Index'); ylabel('Time Sample Index');
saveas(hFigDiffDataContourlet, fullfile(dataFileDir, [dataFileName, '_diffData_contourlet']), 'fig');

fprintf('------------------------------------------------------------\n');
fprintf('Denoised Seismic Data (Contourlet, Thresholding), PSNR = %.2fdB, SSIM = %.4f\n', psnrCleanData_contourlet, ssimCleanData_contourlet);


%% Reference: Contourlet denoising using BPDN with SPG-L1 optimization toolbox
% bpdn:  min ||x||_1  s.t.  ||pdfbFunc(x, 1) - b|| <= \sigma
% lasso: min ||pdfbFunc(x, 1) - b||_2^2 s.t. ||x||_1 < \tau
[vecCoeffContourlet, sContourlet] = pdfb2vec(coeffContourlet);
pdfbFunc = @(x, mode) pdfb(x, sContourlet, pfilter_contourlet, dfilter_contourlet, nlevels_contourlet, nSamples, nRecs, mode);
b = pdfbFunc(vecCoeffContourlet, 1);
% tau = norm(vecCoeffContourlet, 1);
gain = 1;
sigSpThres = sqrt((fp2 - fp1) / (fs/2)) * sigma * sqrt(nSamples * nRecs) * gain;
opts = spgSetParms('verbosity', 1, 'optTol', 1e-6);
% vc_bpdn = spg_lasso(pdfbFunc, b, tau, opts);
vc_bpdn = spg_bpdn(pdfbFunc, b, sigSpThres, opts);

cleanData_contourlet_bpdn = pdfbrec(vec2pdfb(vc_bpdn, sContourlet), pfilter_contourlet, dfilter_contourlet);
save(fullfile(dataFileDir, [dataFileName, '_cleanData_contourlet_bpdn.mat']), 'cleanData_contourlet_bpdn', '-v7.3');
diffData_contourlet_bpdn = dataTrue - cleanData_contourlet_bpdn;
save(fullfile(dataFileDir, [dataFileName, '_diffData_contourlet_bpdn.mat']), 'diffData_contourlet_bpdn', '-v7.3');

% Plot figures and PSNR output
hFigCleanedDataContourlet_bpdn = figure; imagesc(1:nRecs, 1:nSamples, cleanData_contourlet_bpdn); colormap(gray); colorbar; truesize;
psnrCleanData_contourlet_bpdn = 20*log10(sqrt(numel(cleanData_contourlet_bpdn)) / norm(dataTrue(:) - cleanData_contourlet_bpdn(:), 2));
ssimCleanData_contourlet_bpdn = ssim(cleanData_contourlet_bpdn, dataTrue);
xlabel('Trace Index'); ylabel('Time Sample Index');
title(sprintf('Denoised Seismic Data (Contourlet, BPDN), PSNR = %.2fdB', psnrCleanData_contourlet_bpdn));
saveas(hFigCleanedDataContourlet_bpdn, fullfile(dataFileDir, [dataFileName, '_cleanData_contourlet_bpdn']), 'fig');

hFigDiffDataContourlet_bpdn = figure; imagesc(1:nRecs, 1:nSamples, diffData_contourlet_bpdn); colormap(gray); colorbar; truesize;
xlabel('Trace Index'); ylabel('Time Sample Index');
saveas(hFigDiffDataContourlet_bpdn, fullfile(dataFileDir, [dataFileName, '_diffData_contourlet_bpdn']), 'fig');

fprintf('------------------------------------------------------------\n');
fprintf('Denoised Seismic Data (Contourlet, BPDN), PSNR = %.2fdB, SSIM = %.4f\n', psnrCleanData_contourlet_bpdn, ssimCleanData_contourlet_bpdn);


%% Reference: denoising using Curvelet
% 2dB loss if we replace fdct_wrapping and ifdct_wrapping with fdct_usfft
% and ifdct_usfft. respectively
is_real = 1;
nbscales = 4;
nbangles_coarse = 16;

if ~isunix
    coeffCurvelet = fdct_wrapping(noisyData, is_real, 2, nbscales, nbangles_coarse);
else
    coeffCurvelet = fdct_wrapping(noisyData, is_real, nbscales, nbangles_coarse);
end
coeffCurvelet_thres = cell(1, nbscales);

% Set up thresholds for all scales
F = ones(nSamples, nRecs);
X = fftshift(ifft2(F)) * sqrt(nSamples * nRecs);
if ~isunix
    coeffX = fdct_wrapping(X, 0, 2, nbscales, nbangles_coarse);
else
    coeffX = fdct_wrapping(X, 0, nbscales, nbangles_coarse);
end
thres_curvelet = cell(size(coeffX));
for s = 1:length(coeffX)
    thres_curvelet{s} = cell(size(coeffX{s}));
    for w = 1:length(coeffX{s})
        A = coeffX{s}{w};
        thres_curvelet{s}{w} = 3 * sigma * sqrt(sum(sum(A.*conj(A))) / numel(A));
        if s == length(coeffX)
            thres_curvelet{s}{w} = (4/3) * thres_curvelet{s}{w};
        end
    end
end

% Thresholding
for s = 1:length(coeffCurvelet)
    for w = 1:length(coeffCurvelet{s})
        coeffCurvelet_thres{s}{w} = coeffCurvelet{s}{w} .* (abs(coeffCurvelet{s}{w}) > thres_curvelet{s}{w});
    end
end

% Reconstruction
if ~isunix
    cleanData_curvelet = real(ifdct_wrapping(coeffCurvelet_thres, is_real));
else
    cleanData_curvelet = real(ifdct_wrapping(coeffCurvelet_thres, is_real, nbscales, nbangles_coarse));
end
save(fullfile(dataFileDir, [dataFileName, '_cleanData_curvelet.mat']), 'cleanData_curvelet', '-v7.3');
diffData_curvelet = dataTrue - cleanData_curvelet;
save(fullfile(dataFileDir, [dataFileName, '_diffData_curvelet.mat']), 'diffData_curvelet', '-v7.3');

% Plot figures and PSNR output
hFigCleanedDataCurvelet = figure; imagesc(1:nRecs, 1:nSamples, cleanData_curvelet); colormap(gray); colorbar; truesize;
psnrCleanData_curvelet = 20*log10(sqrt(numel(cleanData_curvelet)) / norm(dataTrue(:) - cleanData_curvelet(:), 2));
ssimCleanData_curvelet = ssim(cleanData_curvelet, dataTrue);
xlabel('Trace Index'); ylabel('Time Sample Index');
title(sprintf('Denoised Seismic Data (Curvelet, Thresholding), PSNR = %.2fdB', psnrCleanData_curvelet));
saveas(hFigCleanedDataCurvelet, fullfile(dataFileDir, [dataFileName, '_cleanData_curvelet']), 'fig');

hFigDiffDataCurvelet = figure; imagesc(1:nRecs, 1:nSamples, diffData_curvelet); colormap(gray); colorbar; truesize;
xlabel('Trace Index'); ylabel('Time Sample Index');
saveas(hFigDiffDataCurvelet, fullfile(dataFileDir, [dataFileName, '_diffData_curvelet']), 'fig');

fprintf('------------------------------------------------------------\n');
fprintf('Denoised Seismic Data (Curvelet, Thresholding), PSNR = %.2fdB, SSIM = %.4f\n', psnrCleanData_curvelet, ssimCleanData_curvelet);


%% Reference: Curvelet denoising using BPDN with SPG-L1 optimization toolbox
% bpdn:  min ||x||_1  s.t.  ||fdctFunc(x, 1) - b|| <= \sigma
% lasso: min ||fdctFunc(x, 1) - b||_2^2 s.t. ||x||_1 < \tau
[vecCoeffCurvelet, sCurvelet] = curvelet2vec(coeffCurvelet);
fdctFunc = @(x, mode) fdct(x, sCurvelet, is_real, nbscales, nbangles_coarse, nSamples, nRecs, mode);
b = fdctFunc(vecCoeffCurvelet, 1);
% tau = norm(vecCoeffCurvelet, 1);
gain = 1;
sigSpThres = sqrt((fp2 - fp1) / (fs/2)) * sigma * sqrt(nSamples * nRecs) * gain;
opts = spgSetParms('verbosity', 1, 'optTol', 1e-6);
% vc_bpdn = spg_lasso(fdctFunc, b, tau, opts);
vc_bpdn = spg_bpdn(fdctFunc, b, sigSpThres, opts);

coeffCurvelet_bpdn = vec2curvelet(vc_bpdn, sCurvelet);
if ~isunix
    cleanData_curvelet_bpdn = ifdct_wrapping(coeffCurvelet_bpdn, is_real);
else
    cleanData_curvelet_bpdn = ifdct_wrapping(coeffCurvelet_bpdn, is_real, nbscales, nbangles_coarse);
end
cleanData_curvelet_bpdn = real(cleanData_curvelet_bpdn);
save(fullfile(dataFileDir, [dataFileName, '_cleanData_curvelet_bpdn.mat']), 'cleanData_curvelet_bpdn', '-v7.3');
diffData_curvelet_bpdn = dataTrue - cleanData_curvelet_bpdn;
save(fullfile(dataFileDir, [dataFileName, '_diffData_curvelet_bpdn.mat']), 'diffData_curvelet_bpdn', '-v7.3');

% Plot figures and PSNR output
hFigCleanedDataCurvelet_bpdn = figure; imagesc(1:nRecs, 1:nSamples, cleanData_curvelet_bpdn); colormap(gray); colorbar; truesize;
psnrCleanData_curvelet_bpdn = 20*log10(sqrt(numel(cleanData_curvelet_bpdn)) / norm(dataTrue(:) - cleanData_curvelet_bpdn(:), 2));
ssimCleanData_curvelet_bpdn = ssim(cleanData_curvelet_bpdn, dataTrue);
xlabel('Trace Index'); ylabel('Time Sample Index');
title(sprintf('Denoised Seismic Data (Curvelet, BPDN), PSNR = %.2fdB', psnrCleanData_curvelet_bpdn));
saveas(hFigCleanedDataCurvelet_bpdn, fullfile(dataFileDir, [dataFileName, '_cleanData_curvelet_bpdn']), 'fig');

hFigDiffDataCurvelet_bpdn = figure; imagesc(1:nRecs, 1:nSamples, diffData_curvelet_bpdn); colormap(gray); colorbar; truesize;
xlabel('Trace Index'); ylabel('Time Sample Index');
saveas(hFigDiffDataCurvelet_bpdn, fullfile(dataFileDir, [dataFileName, '_diffData_curvelet_bpdn']), 'fig');

fprintf('------------------------------------------------------------\n');
fprintf('Denoised Seismic Data (Curvelet, BPDN), PSNR = %.2fdB, SSIM = %.4f\n', psnrCleanData_curvelet_bpdn, ssimCleanData_curvelet_bpdn);


%% %% Dictionary learning using sparse K-SVD
close all;
%% Base dictionary setting
nlevels_wavelet_train = [0, 0];       % Decomposition level, all 0 means wavelet
pfilter_wavelet_train = '9/7' ;       % Pyramidal filter
dfilter_wavelet_train = 'pkva' ;      % Directional filter

is_real_train = 1;
nbscales_train = 3;
nbangles_coarse_train = 16;


% %% Dictionary learning using sparse K-SVD (band-by-band)
% % wavelet transform of train data
% trainData_coeff = pdfbdec(trainData, pfilter_wavelet_train, dfilter_wavelet_train, nlevels_wavelet_train);
% cleanCoeff_sparseKsvd = cell(size(trainData_coeff));
% learnedDict_band = cell(size(trainData_coeff));
% learnedDictImg_band = cell(size(trainData_coeff));
% overallDict_band = cell(size(trainData_coeff));
% overallDictImg_band = cell(size(trainData_coeff));
% nLayer = length(trainData_coeff);
% for iLayer = 2:nLayer
%     cleanCoeff_sparseKsvd{iLayer} = cell(size(trainData_coeff{iLayer}));
%     learnedDict_band{iLayer} = cell(size(trainData_coeff{iLayer}));
%     learnedDictImg_band{iLayer} = cell(size(trainData_coeff{iLayer}));
%     overallDict_band{iLayer} = cell(size(trainData_coeff{iLayer}));
%     overallDictImg_band{iLayer} = cell(size(trainData_coeff{iLayer}));
% end
% fprintf('------------------------------------------------------------\n');
% fprintf('Dictionary Learning (band-by-band)\n');
% idxBand = 0;
% % coarest band
% % parameters
% gain_band = 1;                                  % noise gain (default value 1.15)
% trainBlockSize_band = 16;                        % for each dimension
% trainBlockNum_band = 5000;                      % number of training blocks in the training set
% trainIter_band = 20;
% sigSpThres_band = sqrt((fp2 - fp1) / (fs/2)) * sigma * trainBlockSize_band * gain_band / 2^length(nlevels_wavelet_train); % pre-defined l2-norm error for BPDN
% atomSpThres_band = 50;                          % a self-determind value to control the sparsity of matrix A
% % dictionary learning
% fprintf('Learning from Band %d\n', idxBand);
% initDict_band = speye(trainBlockSize_band^2, trainBlockSize_band^2);
% % test begin
% % wavelets for base dictionary
% [vecTrainBlockCoeff_band, str_band] = pdfb2vec(pdfbdec(zeros(trainBlockSize_band, trainBlockSize_band), pfilter_wavelet_train, dfilter_wavelet_train, nlevels_wavelet_train));
% baseSynOp_band = @(x) pdfb(x, str_band, pfilter_wavelet_train, dfilter_wavelet_train, nlevels_wavelet_train, trainBlockSize_band, trainBlockSize_band, 1);
% baseAnaOp_band = @(x) pdfb(x, str_band, pfilter_wavelet_train, dfilter_wavelet_train, nlevels_wavelet_train, trainBlockSize_band, trainBlockSize_band, 2);
% [PhiSyn_band, PhiAna_band] = operator2matrix(baseSynOp_band, baseAnaOp_band, trainBlockSize_band * trainBlockSize_band);
% PhiSyn_band = PhiSyn_band(:, 1 : prod(str_band(1, 3:4)));
% PhiAna_band = PhiAna_band(1 : prod(str_band(1, 3:4)), :);
% idxAtom = prod(str_band(1, 3:4));
% initDict_band = speye(prod(str_band(1, 3:4)), prod(str_band(1, 3:4)));
% % test end
% [learnedDict_band{1}, Coeffs_band, errLasso_band, errBpdn_band] = sparseKsvd(trainData_coeff{1}, PhiSyn_band, PhiAna_band, ...
%     initDict_band, trainIter_band, trainBlockSize_band, trainBlockNum_band, atomSpThres_band, sigSpThres_band, struct('verbosity', 0, 'method', 'bpdn'));
% % test begin
% learnedDictImg_band{1} = showdict(full(learnedDict_band{1}), [1 1]*sqrt(size(learnedDict_band{1}, 1)), round(sqrt(size(learnedDict_band{1}, 2))), round(sqrt(size(learnedDict_band{1}, 2))), 'whitelines', 'highcontrast');
% overallDict_band{1} = PhiSyn_band * learnedDict_band{1};
% overallDictImg_band{1} = showdict(PhiSyn_band * learnedDict_band{1}, [1 1]*sqrt(size(PhiSyn_band * learnedDict_band{1}, 1)), round(sqrt(size(PhiSyn_band * learnedDict_band{1}, 2))), round(sqrt(size(PhiSyn_band * learnedDict_band{1}, 2))), 'whitelines', 'highcontrast');
% % test end
% % % start a pool of Matlab workers
% % numCores = feature('numcores');
% % if isempty(gcp('nocreate')) % checking to see if my pool is already open
% %     myPool = parpool(numCores);
% % end
% % % denoising
% % fprintf('Denoising Band %d\n', idxBand);
% % cleanCoeff_sparseKsvd{1} = zeros(size(trainData_coeff{1}));
% % totalBlockNum_band = (size(trainData_coeff{1}, 1) - trainBlockSize_band + 1) * (size(trainData_coeff{1}, 2) - trainBlockSize_band + 1);
% % processedBlocks = 0;
% % for ibatch = 1:size(trainData_coeff{1}, 2)-trainBlockSize_band+1
% %     fprintf('\tBatch %d... ', ibatch);
% %     % the current batch of blocks
% %     blocks = im2colstep(trainData_coeff{1}(:, ibatch:ibatch+trainBlockSize_band-1), trainBlockSize_band * [1, 1], [1, 1]);
% %     cleanBlocks = zeros(size(blocks));
% %     blockCoeff = zeros(trainBlockSize_band^2, size(trainData_coeff{1}, 1) - trainBlockSize_band + 1);
% %     parfor iblk = 1:size(trainData_coeff{1}, 1) - trainBlockSize_band + 1
% %         opts = spgSetParms('verbosity', 0, 'optTol', 1e-6);
% %         blockCoeff(:, iblk) = spg_bpdn(@(x, mode) learnedOp(x, baseSynOp_band, baseAnaOp_band, learnedDict_band, mode), blocks(:, iblk), sigSpThres_band, opts);
% %         % blockCoeff(:, iblk) = OMP({@(x) baseSynOp(learnedDict*x), @(x) learnedDict'*baseAnaOp(x)}, blocks(:, iblk), sigSpThres);
% %         cleanBlocks(:, iblk) = learnedOp(blockCoeff(:, iblk), baseSynOp_band, baseAnaOp_band, learnedDict_band, 1);
% %     end
% %     cleanBatch = col2imstep(cleanBlocks, [size(trainData_coeff{1}, 1), trainBlockSize_band], trainBlockSize_band * [1, 1], [1, 1]);
% %     cleanCoeff_sparseKsvd{1}(:,ibatch:ibatch+trainBlockSize_band-1) = cleanCoeff_sparseKsvd{1}(:,ibatch:ibatch+trainBlockSize_band-1) + cleanBatch;
% %     processedBlocks = processedBlocks + (size(trainData_coeff{1}, 1) - trainBlockSize_band + 1);
% %     fprintf('\tProcessed %d blocks\n', processedBlocks);
% % end
% % % average the denoised and noisy signals
% % cnt = countcover(size(trainData_coeff{1}), trainBlockSize_band * [1, 1], [1, 1]);
% % cleanCoeff_sparseKsvd{1} = cleanCoeff_sparseKsvd{1}./cnt;
% % % terminate the pool of Matlab workers
% % delete(gcp('nocreate'));
% % other bands
% for iLayer = 2:nLayer
%     nDir = length(trainData_coeff{iLayer});
%     for iDir = 1:nDir
%         idxBand = idxBand + 1;
%         % dictionary learning
%         fprintf('Learning from Band %d\n', idxBand);
%         initDict_band = speye(trainBlockSize_band^2, trainBlockSize_band^2);
%         % test begin
%         % wavelets for base dictionary
%         [vecTrainBlockCoeff_band, str_band] = pdfb2vec(pdfbdec(zeros(trainBlockSize_band, trainBlockSize_band), pfilter_wavelet_train, dfilter_wavelet_train, nlevels_wavelet_train));
%         baseSynOp_band = @(x) pdfb(x, str_band, pfilter_wavelet_train, dfilter_wavelet_train, nlevels_wavelet_train, trainBlockSize_band, trainBlockSize_band, 1);
%         baseAnaOp_band = @(x) pdfb(x, str_band, pfilter_wavelet_train, dfilter_wavelet_train, nlevels_wavelet_train, trainBlockSize_band, trainBlockSize_band, 2);
%         [PhiSyn_band, PhiAna_band] = operator2matrix(baseSynOp_band, baseAnaOp_band, trainBlockSize_band * trainBlockSize_band);
%         PhiSyn_band = PhiSyn_band(:, idxAtom + 1 : idxAtom + prod(str_band(1+(iLayer-2)*nDir+iDir, 3:4)));
%         PhiAna_band = PhiAna_band(idxAtom + 1 : idxAtom + prod(str_band(1+(iLayer-2)*nDir+iDir, 3:4)), :);
%         idxAtom = idxAtom + prod(str_band(1+(iLayer-2)*nDir+iDir, 3:4));
%         initDict_band = speye(prod(str_band(1+(iLayer-2)*nDir+iDir, 3:4)), prod(str_band(1+(iLayer-2)*nDir+iDir, 3:4)));
%         % test end
%         [learnedDict_band{iLayer}{iDir}, Coeffs_band, errLasso_band, errBpdn_band] = sparseKsvd(trainData_coeff{iLayer}{iDir}, PhiSyn_band, PhiAna_band, ...
%             initDict_band, trainIter_band, trainBlockSize_band, trainBlockNum_band, atomSpThres_band, sigSpThres_band, struct('verbosity', 0, 'method', 'bpdn'));
%         % test begin
%         learnedDictImg_band{iLayer}{iDir} = showdict(full(learnedDict_band{iLayer}{iDir}), [1 1]*sqrt(size(learnedDict_band{iLayer}{iDir}, 1)), round(sqrt(size(learnedDict_band{iLayer}{iDir}, 2))), round(sqrt(size(learnedDict_band{iLayer}{iDir}, 2))), 'whitelines', 'highcontrast');
%         overallDict_band{iLayer}{iDir} = PhiSyn_band * learnedDict_band{iLayer}{iDir};
%         overallDictImg_band{iLayer}{iDir} = showdict(PhiSyn_band * learnedDict_band{iLayer}{iDir}, [1 1]*sqrt(size(PhiSyn_band * learnedDict_band{iLayer}{iDir}, 1)), round(sqrt(size(PhiSyn_band * learnedDict_band{iLayer}{iDir}, 2))), round(sqrt(size(PhiSyn_band * learnedDict_band{iLayer}{iDir}, 2))), 'whitelines', 'highcontrast');
%         % test end
% %         % start a pool of Matlab workers
% %         numCores = feature('numcores');
% %         if isempty(gcp('nocreate')) % checking to see if my pool is already open
% %             myPool = parpool(numCores);
% %         end
% %         % denoising
% %         fprintf('Denoising Band %d\n', idxBand);
% %         cleanCoeff_sparseKsvd{iLayer}{iDir} = zeros(size(trainData_coeff{iLayer}{iDir}));
% %         totalBlockNum_band = (size(trainData_coeff{iLayer}{iDir}, 1) - trainBlockSize_band + 1) * (size(trainData_coeff{iLayer}{iDir}, 2) - trainBlockSize_band + 1);
% %         processedBlocks = 0;
% %         for ibatch = 1:size(trainData_coeff{iLayer}{iDir}, 2)-trainBlockSize_band+1
% %             fprintf('\tBatch %d... ', ibatch);
% %             % the current batch of blocks
% %             blocks = im2colstep(trainData_coeff{iLayer}{iDir}(:, ibatch:ibatch+trainBlockSize_band-1), trainBlockSize_band * [1, 1], [1, 1]);
% %             cleanBlocks = zeros(size(blocks));
% %             blockCoeff = zeros(trainBlockSize_band^2, size(trainData_coeff{iLayer}{iDir}, 1) - trainBlockSize_band + 1);
% %             parfor iblk = 1:size(trainData_coeff{iLayer}{iDir}, 1) - trainBlockSize_band + 1
% %                 opts = spgSetParms('verbosity', 0, 'optTol', 1e-6);
% %                 blockCoeff(:, iblk) = spg_bpdn(@(x, mode) learnedOp(x, baseSynOp_band, baseAnaOp_band, learnedDict_band, mode), blocks(:, iblk), sigSpThres_band, opts);
% %                 % blockCoeff(:, iblk) = OMP({@(x) baseSynOp(learnedDict*x), @(x) learnedDict'*baseAnaOp(x)}, blocks(:, iblk), sigSpThres);
% %                 cleanBlocks(:, iblk) = learnedOp(blockCoeff(:, iblk), baseSynOp_band, baseAnaOp_band, learnedDict_band, 1);
% %             end
% %             cleanBatch = col2imstep(cleanBlocks, [size(trainData_coeff{iLayer}{iDir}, 1), trainBlockSize_band], trainBlockSize_band * [1, 1], [1, 1]);
% %             cleanCoeff_sparseKsvd{iLayer}{iDir}(:,ibatch:ibatch+trainBlockSize_band-1) = cleanCoeff_sparseKsvd{iLayer}{iDir}(:,ibatch:ibatch+trainBlockSize_band-1) + cleanBatch;
% %             processedBlocks = processedBlocks + (size(trainData_coeff{iLayer}{iDir}, 1) - trainBlockSize_band + 1);
% %             fprintf('\tProcessed %d blocks\n', processedBlocks);
% %         end
% %         % average the denoised and noisy signals
% %         cnt = countcover(size(trainData_coeff{iLayer}{iDir}), trainBlockSize_band * [1, 1], [1, 1]);
% %         cleanCoeff_sparseKsvd{iLayer}{iDir} = cleanCoeff_sparseKsvd{iLayer}{iDir}./cnt;
% %         % terminate the pool of Matlab workers
% %         delete(gcp('nocreate'));
%     end
%     sigSpThres_band = sigSpThres_band * 2;
% end
% % test begin
% save(fullfile(dataFileDir, [dataFileName, '_learnedDict_band.mat']), 'learnedDict_band', '-v7.3');
% save(fullfile(dataFileDir, [dataFileName, '_learnedDictImg_band.mat']), 'learnedDictImg_band', '-v7.3');
% save(fullfile(dataFileDir, [dataFileName, '_overallDict_band.mat']), 'overallDict_band', '-v7.3');
% save(fullfile(dataFileDir, [dataFileName, '_overallDictImg_band.mat']), 'overallDictImg_band', '-v7.3');
% % test end
% % % reconstruction
% % cleanData_sparseKsvd_band = pdfbrec(cleanCoeff_sparseKsvd, pfilter_wavelet_train, dfilter_wavelet_train);
% % save(fullfile(dataFileDir, [dataFileName, '_cleanData_sparseKsvd_band.mat']), 'cleanData_sparseKsvd_band', '-v7.3');
% % diffData_sparseKsvd_band = dataTrue - cleanData_sparseKsvd_band;
% % save(fullfile(dataFileDir, [dataFileName, '_diffData_sparseKsvd_band.mat']), 'diffData_sparseKsvd_band', '-v7.3');
% 
% % %% Plot figures and PSNR output
% % hFigCleanedDataSparseKsvd_band = figure; imagesc(1:nRecs, 1:nSamples, cleanData_sparseKsvd_band); colormap(gray); colorbar; truesize;
% % psnrCleanData_sparseKsvd_band = 20*log10(sqrt(numel(cleanData_sparseKsvd_band)) / norm(dataTrue(:) - cleanData_sparseKsvd_band(:), 2));
% % ssimCleanData_sparseKsvd_band = ssim(cleanData_sparseKsvd_band, dataTrue);
% % xlabel('Trace Index'); ylabel('Time Sample Index');
% % title(sprintf('Denoised Seismic Data (Dictionary Learning, band-by-band), PSNR = %.2fdB', psnrCleanData_sparseKsvd_band));
% % saveas(hFigCleanedDataSparseKsvd_band, fullfile(dataFileDir, [dataFileName, '_cleanData_sparseKsvd_band']), 'fig');
% % 
% % hFigDiffDataSparseKsvd_band = figure; imagesc(1:nRecs, 1:nSamples, diffData_sparseKsvd_band); colormap(gray); colorbar; truesize;
% % xlabel('Trace Index'); ylabel('Time Sample Index');
% % saveas(hFigDiffDataSparseKsvd_band, fullfile(dataFileDir, [dataFileName, '_diffData_sparseKsvd_band']), 'fig');
% % 
% % fprintf('------------------------------------------------------------\n');
% % fprintf('Denoised Seismic Data (Dictionary Learning, band-by-band), PSNR = %.2fdB, SSIM = %.4f\n', psnrCleanData_sparseKsvd_band, ssimCleanData_sparseKsvd_band);
% 

%% Dictionary learning using sparse K-SVD (all bands together)
fprintf('------------------------------------------------------------\n');
fprintf('Dictionary Learning (all bands together)\n');
gain = 1;                                   % noise gain (default value 1.15)
trainBlockSize = 16;                        % for each dimension
trainBlockNum = 5000;                       % number of training blocks in the training set
trainIter = 20;
sigSpThres = sqrt((fp2 - fp1) / (fs/2)) * sigma * trainBlockSize * gain; % pre-defined l2-norm error for BPDN
atomSpThres = 50;                          % a self-determind value to control the sparsity of matrix A

% wavelets for base dictionary
[vecTrainBlockCoeff, str] = pdfb2vec(pdfbdec(zeros(trainBlockSize, trainBlockSize), pfilter_wavelet_train, dfilter_wavelet_train, nlevels_wavelet_train));
initDict = speye(length(vecTrainBlockCoeff), length(vecTrainBlockCoeff));
baseSynOp = @(x) pdfb(x, str, pfilter_wavelet_train, dfilter_wavelet_train, nlevels_wavelet_train, trainBlockSize, trainBlockSize, 1);
baseAnaOp = @(x) pdfb(x, str, pfilter_wavelet_train, dfilter_wavelet_train, nlevels_wavelet_train, trainBlockSize, trainBlockSize, 2);

% % curvelets for base dictionary
% if ~isunix
%     [vecTrainBlockCoeff, str] = curvelet2vec(fdct_wrapping(zeros(trainBlockSize, trainBlockSize), is_real_train, 2, nbscales_train, nbangles_coarse_train));
% else
%     [vecTrainBlockCoeff, str] = curvelet2vec(fdct_wrapping(zeros(trainBlockSize, trainBlockSize), is_real_train, nbscales_train, nbangles_coarse_train));
% end
% initDict = eye(length(vecTrainBlockCoeff), length(vecTrainBlockCoeff));
% baseSynOp = @(x) fdct(x, str, is_real_train, nbscales_train, nbangles_coarse_train, trainBlockSize, trainBlockSize, 1);
% baseAnaOp = @(x) fdct(x, str, is_real_train, nbscales_train, nbangles_coarse_train, trainBlockSize, trainBlockSize, 2);

[learnedDict, Coeffs, errLasso, errBpdn] = sparseKsvd(trainData, baseSynOp, baseAnaOp, ...
    initDict, trainIter, trainBlockSize, trainBlockNum, atomSpThres, sigSpThres, struct('verbosity', 0, 'method', 'bpdn'));


%% Show trained dictionary
if (isa(baseSynOp, 'function_handle') && isa(baseAnaOp, 'function_handle'))
    [PhiSyn, PhiAna] = operator2matrix(baseSynOp, baseAnaOp, trainBlockSize * trainBlockSize);
else
    PhiSyn = baseSynOp;
    PhiAna = baseAnaOp;
end
% show atoms of base dictionary
PhiSynImg = showdict(PhiSyn, [1 1]*sqrt(size(PhiSyn, 1)), round(sqrt(size(PhiSyn, 2))), round(sqrt(size(PhiSyn, 2))), 'whitelines', 'highcontrast');
hFigPhiSyn = figure; imshow(imresize(PhiSynImg, 2, 'nearest')); title('Daubechies 9/7 wavelets');
saveas(hFigPhiSyn, fullfile(dataFileDir, [dataFileName, '_db97_wavelet']), 'fig');
% show atoms A
AImg = showdict(full(learnedDict), [1 1]*sqrt(size(learnedDict, 1)), round(sqrt(size(learnedDict, 2))), round(sqrt(size(learnedDict, 2))), 'whitelines', 'highcontrast');
hFigA = figure; imshow(imresize(AImg, 2, 'nearest')); title(sprintf('Atoms in Matrix A (%d iterations)', trainIter));
saveas(hFigA, fullfile(dataFileDir, [dataFileName, '_A']), 'fig');
% show atoms of overall learned dictionary
dictImg = showdict(PhiSyn * learnedDict, [1 1]*sqrt(size(PhiSyn * learnedDict, 1)), round(sqrt(size(PhiSyn * learnedDict, 2))), round(sqrt(size(PhiSyn * learnedDict, 2))), 'whitelines', 'highcontrast');
hFigLearnedDict = figure; imshow(imresize(dictImg, 2, 'nearest')); title(sprintf('Trained Dictionary (%d iterations)', trainIter));
saveas(hFigLearnedDict, fullfile(dataFileDir, [dataFileName, '_learnedDict']), 'fig');


%% Start a pool of Matlab workers
numCores = feature('numcores');
if isempty(gcp('nocreate')) % checking to see if my pool is already open
    myPool = parpool(numCores);
end


%% Denoising
fprintf('------------------------------------------------------------\n');
fprintf('Denoising\n');

cleanData_sparseKsvd = zeros(size(dataTrue));
totalBlockNum = (nSamples - trainBlockSize + 1) * (nRecs - trainBlockSize + 1);
processedBlocks = 0;

for ibatch = 1:nRecs-trainBlockSize+1
    fprintf('Batch %d... ', ibatch);
    % the current batch of blocks
    blocks = im2colstep(noisyData(:, ibatch:ibatch+trainBlockSize-1), trainBlockSize * [1, 1], [1, 1]);
    
    % % remove DC (mean values)
    % [blocks, dc] = remove_dc(blocks,'columns');
    
    cleanBlocks = zeros(size(blocks));
    blockCoeff = zeros(length(vecTrainBlockCoeff), nSamples - trainBlockSize + 1);
    parfor iblk = 1:nSamples - trainBlockSize + 1
        opts = spgSetParms('verbosity', 0, 'optTol', 1e-6);
        blockCoeff(:, iblk) = spg_bpdn(@(x, mode) learnedOp(x, baseSynOp, baseAnaOp, learnedDict, mode), blocks(:, iblk), sigSpThres, opts);
        % blockCoeff(:, iblk) = OMP({@(x) baseSynOp(learnedDict*x), @(x) learnedDict'*baseAnaOp(x)}, blocks(:, iblk), sigSpThres);
        cleanBlocks(:, iblk) = learnedOp(blockCoeff(:, iblk), baseSynOp, baseAnaOp, learnedDict, 1);
    end
    
    % % add DC (mean values)
    % cleanBlocks = add_dc(cleanBlocks, dc, 'columns');
    
    cleanBatch = col2imstep(cleanBlocks, [nSamples, trainBlockSize], trainBlockSize * [1, 1], [1, 1]);
    cleanData_sparseKsvd(:,ibatch:ibatch+trainBlockSize-1) = cleanData_sparseKsvd(:,ibatch:ibatch+trainBlockSize-1) + cleanBatch;
    
    processedBlocks = processedBlocks + (nSamples - trainBlockSize + 1);
    fprintf('Processed %d blocks\n', processedBlocks);
end

% average the denoised and noisy signals
cnt = countcover(size(trainData), trainBlockSize * [1, 1], [1, 1]);
cleanData_sparseKsvd = cleanData_sparseKsvd./cnt;
save(fullfile(dataFileDir, [dataFileName, '_cleanData_sparseKsvd.mat']), 'cleanData_sparseKsvd', '-v7.3');
diffData_sparseKsvd = dataTrue - cleanData_sparseKsvd;
save(fullfile(dataFileDir, [dataFileName, '_diffData_sparseKsvd.mat']), 'diffData_sparseKsvd', '-v7.3');

%% Terminate the pool of Matlab workers
delete(gcp('nocreate'));

%% Plot figures and PSNR output
hFigCleanedDataSparseKsvd = figure; imagesc(1:nRecs, 1:nSamples, cleanData_sparseKsvd); colormap(gray); colorbar; truesize;
psnrCleanData_sparseKsvd = 20*log10(sqrt(numel(cleanData_sparseKsvd)) / norm(dataTrue(:) - cleanData_sparseKsvd(:), 2));
ssimCleanData_sparseKsvd = ssim(cleanData_sparseKsvd, dataTrue);
xlabel('Trace Index'); ylabel('Time Sample Index');
title(sprintf('Denoised Seismic Data (Dictionary Learning, all bands together), PSNR = %.2fdB', psnrCleanData_sparseKsvd));
saveas(hFigCleanedDataSparseKsvd, fullfile(dataFileDir, [dataFileName, '_cleanData_sparseKsvd']), 'fig');

hFigDiffDataSparseKsvd = figure; imagesc(1:nRecs, 1:nSamples, diffData_sparseKsvd); colormap(gray); colorbar; truesize;
xlabel('Trace Index'); ylabel('Time Sample Index');
saveas(hFigDiffDataSparseKsvd, fullfile(dataFileDir, [dataFileName, '_diffData_sparseKsvd']), 'fig');

fprintf('------------------------------------------------------------\n');
fprintf('Denoised Seismic Data (Dictionary Learning, all bands together), PSNR = %.2fdB, SSIM = %.4f\n', psnrCleanData_sparseKsvd, ssimCleanData_sparseKsvd);
