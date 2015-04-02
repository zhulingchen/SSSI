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
xlabel('Traces'); ylabel('Time');
title('Original Seismic Data');

hFigNoisyData = figure; imagesc(1:nRecs, 1:nSamples, noisyData); colormap(gray); colorbar; truesize;
psnrNoisyData = 20*log10(sqrt(numel(noisyData)) / norm(dataTrue(:) - noisyData(:), 2));
ssimNoisyData = ssim(noisyData, dataTrue);
xlabel('Traces'); ylabel('Time');
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
xlabel('Traces'); ylabel('Time');
title(sprintf('Denoised Seismic Data (Wavelet), PSNR = %.2fdB', psnrCleanData_wavelet));
saveas(hFigCleanedDataWavelet, fullfile(dataFileDir, [dataFileName, '_cleanData_wavelet']), 'fig');

hFigDiffDataWavelet = figure; imagesc(1:nRecs, 1:nSamples, diffData_wavelet); colormap(gray); colorbar; truesize;
xlabel('Traces'); ylabel('Time');
saveas(hFigDiffDataWavelet, fullfile(dataFileDir, [dataFileName, '_diffData_wavelet']), 'fig');

fprintf('------------------------------------------------------------\n');
fprintf('Denoised Seismic Data (Wavelet), PSNR = %.2fdB, SSIM = %.4f\n', psnrCleanData_wavelet, ssimCleanData_wavelet);


%% Reference: denoising using Contourlet
nlevels_contourlet = [2, 3, 4];     % Decomposition level, all 0 means wavelet
pfilter_contourlet = '9/7';         % Pyramidal filter
dfilter_contourlet = 'pkva';        % Directional filter
[vecContourletCoeff, str] = pdfb2vec(pdfbdec(noisyData, pfilter_contourlet, dfilter_contourlet, nlevels_contourlet));

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
xlabel('Traces'); ylabel('Time');
title(sprintf('Denoised Seismic Data (Contourlet), PSNR = %.2fdB', psnrCleanData_contourlet));
saveas(hFigCleanedDataContourlet, fullfile(dataFileDir, [dataFileName, '_cleanData_contourlet']), 'fig');

hFigDiffDataContourlet = figure; imagesc(1:nRecs, 1:nSamples, diffData_contourlet); colormap(gray); colorbar; truesize;
xlabel('Traces'); ylabel('Time');
saveas(hFigDiffDataContourlet, fullfile(dataFileDir, [dataFileName, '_diffData_contourlet']), 'fig');

fprintf('------------------------------------------------------------\n');
fprintf('Denoised Seismic Data (Contourlet), PSNR = %.2fdB, SSIM = %.4f\n', psnrCleanData_contourlet, ssimCleanData_contourlet);


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
xlabel('Traces'); ylabel('Time');
title(sprintf('Denoised Seismic Data (Curvelet), PSNR = %.2fdB', psnrCleanData_curvelet));
saveas(hFigCleanedDataCurvelet, fullfile(dataFileDir, [dataFileName, '_cleanData_curvelet']), 'fig');

hFigDiffDataCurvelet = figure; imagesc(1:nRecs, 1:nSamples, diffData_curvelet); colormap(gray); colorbar; truesize;
xlabel('Traces'); ylabel('Time');
saveas(hFigDiffDataCurvelet, fullfile(dataFileDir, [dataFileName, '_diffData_curvelet']), 'fig');

fprintf('------------------------------------------------------------\n');
fprintf('Denoised Seismic Data (Curvelet), PSNR = %.2fdB, SSIM = %.4f\n', psnrCleanData_curvelet, ssimCleanData_curvelet);


% Curvelet denoising using BPDN with SPG-L1 optimization toolbox
% close all;
% bpdn:  min ||x||_1  s.t.  ||fdctFunc(x, 1) - b|| <= \sigma
% lasso: min ||fdctFunc(x, 1) - b||_2^2 s.t. ||x||_1 < \tau
[vecCoeffCurvelet, sCurvelet] = curvelet2vec(coeffCurvelet);
fdctFunc = @(x, mode) fdct(x, sCurvelet, is_real, nbscales, nbangles_coarse, nSamples, nRecs, mode);
b = fdctFunc(vecCoeffCurvelet, 1);
% tau = norm(vecCoeffCurvelet, 1);
gain = 1;
sigSpThres = sqrt((fp2 - fp1) / (fs/2)) * sigma * sqrt(nSamples * nRecs) * gain;
opts = spgSetParms('verbosity', 1, 'optTol', 1e-6);
% vc_spg = spg_lasso(fdctFunc, b, tau, opts);
vc_spg = spg_bpdn(fdctFunc, b, sigSpThres, opts);

coeffCurvelet_spg = vec2curvelet(vc_spg, sCurvelet);
if ~isunix
    cleanData_curvelet_spg = ifdct_wrapping(coeffCurvelet_spg, is_real);
else
    cleanData_curvelet_spg = ifdct_wrapping(coeffCurvelet_spg, is_real, nbscales, nbangles_coarse);
end
cleanData_curvelet_spg = real(cleanData_curvelet_spg);
diffData_curvelet_spg = dataTrue - cleanData_curvelet_spg;
save(fullfile(dataFileDir, [dataFileName, '_diffData_curvelet_spg.mat']), 'diffData_curvelet_spg', '-v7.3');

% Plot figures and PSNR output
hFigCleanedDataCurvelet_spg = figure; imagesc(1:nRecs, 1:nSamples, cleanData_curvelet_spg); colormap(gray); colorbar; truesize;
psnrCleanData_curvelet_spg = 20*log10(sqrt(numel(cleanData_curvelet_spg)) / norm(dataTrue(:) - cleanData_curvelet_spg(:), 2));
ssimCleanData_curvelet_spg = ssim(cleanData_curvelet_spg, dataTrue);
xlabel('Traces'); ylabel('Time');
title(sprintf('Denoised Seismic Data (Curvelet, BPDN), PSNR = %.2fdB', psnrCleanData_curvelet_spg));
saveas(hFigCleanedDataCurvelet_spg, fullfile(dataFileDir, [dataFileName, '_cleanData_curvelet_spg']), 'fig');

hFigDiffDataCurvelet_spg = figure; imagesc(1:nRecs, 1:nSamples, diffData_curvelet_spg); colormap(gray); colorbar; truesize;
xlabel('Traces'); ylabel('Time');
saveas(hFigDiffDataCurvelet_spg, fullfile(dataFileDir, [dataFileName, '_diffData_curvelet_spg']), 'fig');

fprintf('------------------------------------------------------------\n');
fprintf('Denoised Seismic Data (Curvelet, BPDN), PSNR = %.2fdB, SSIM = %.4f\n', psnrCleanData_curvelet_spg, ssimCleanData_curvelet_spg);


%% Parameters for dictionary learning using sparse K-SVD
gain = 1;                                   % noise gain (default value 1.15)
trainBlockSize = 16;                        % for each dimension
trainBlockNum = 5000;                       % number of training blocks in the training set
trainIter = 20;
sigSpThres = sqrt((fp2 - fp1) / (fs/2)) * sigma * trainBlockSize * gain; % pre-defined l2-norm error for BPDN
atomSpThres = 50;                          % a self-determind value to control the sparsity of matrix A


%% Base dictionary setting
nlevels_wavelet_train = [0, 0];       % Decomposition level, all 0 means wavelet
pfilter_wavelet_train = '9/7' ;       % Pyramidal filter
dfilter_wavelet_train = 'pkva' ;      % Directional filter

is_real_train = 1;
nbscales_train = 3;
nbangles_coarse_train = 16;


%% Dictionary learning using sparse K-SVD
fprintf('------------------------------------------------------------\n');
fprintf('Dictionary Learning\n');

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

% % wavelet transform of train data
% dataTrue_coeff = pdfbdec(dataTrue, pfilter_wavelet_train, dfilter_wavelet_train, nlevels_wavelet_train);


[learnedDict, Coeffs, errLasso, errBpdn] = sparseKsvd(trainData, baseSynOp, baseAnaOp, ...
    initDict, trainIter, trainBlockSize, trainBlockNum, atomSpThres, sigSpThres, struct('verbosity', 0, 'method', 'bpdn'));


%% show trained dictionary
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
cnt = countcover(size(noisyData), trainBlockSize * [1, 1], [1, 1]);
cleanData_sparseKsvd = cleanData_sparseKsvd./cnt;
save(fullfile(dataFileDir, [dataFileName, '_cleanData_sparseKsvd.mat']), 'cleanData_sparseKsvd', '-v7.3');
diffData_sparseKsvd = dataTrue - cleanData_sparseKsvd;
save(fullfile(dataFileDir, [dataFileName, '_diffData_sparseKsvd.mat']), 'diffData_sparseKsvd', '-v7.3');


%% Plot figures and PSNR output
hFigCleanedDataSparseKsvd = figure; imagesc(1:nRecs, 1:nSamples, cleanData_sparseKsvd); colormap(gray); colorbar; truesize;
psnrCleanData_sparseKsvd = 20*log10(sqrt(numel(cleanData_sparseKsvd)) / norm(dataTrue(:) - cleanData_sparseKsvd(:), 2));
ssimCleanData_sparseKsvd = ssim(cleanData_sparseKsvd, dataTrue);
xlabel('Traces'); ylabel('Time');
title(sprintf('Denoised Seismic Data, PSNR = %.2fdB', psnrCleanData_sparseKsvd));
saveas(hFigCleanedDataSparseKsvd, fullfile(dataFileDir, [dataFileName, '_cleanData_sparseKsvd']), 'fig');

hFigDiffDataSparseKsvd = figure; imagesc(1:nRecs, 1:nSamples, diffData_sparseKsvd); colormap(gray); colorbar; truesize;
xlabel('Traces'); ylabel('Time');
saveas(hFigDiffDataSparseKsvd, fullfile(dataFileDir, [dataFileName, '_diffData_sparseKsvd']), 'fig');

fprintf('------------------------------------------------------------\n');
fprintf('Denoised Seismic Data, PSNR = %.2fdB, SSIM = %.4f\n', psnrCleanData_sparseKsvd, ssimCleanData_sparseKsvd);


%% Terminate the pool of Matlab workers
delete(gcp('nocreate'));
