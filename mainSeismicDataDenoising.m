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

load('./modelData/denoising/dataTrue.mat'); % dataTrue
nBoundary = 20;
dataTrue = dataTrue(nBoundary+1:end-nBoundary,:)';
[nSamples, nRecs] = size(dataTrue);

%% Remove mean and normalization
minData = min(dataTrue, [], 1);
maxData = max(dataTrue, [], 1);
meanData = mean(dataTrue, 1);
% dataTrue = bsxfun(@times, bsxfun(@minus, dataTrue, minData), 1./abs(maxData - minData));
dataTrue = bsxfun(@times, bsxfun(@minus, dataTrue, meanData), 1./abs(maxData - minData));

% %% Normalize data to unit norm
% for ir = 1:nRecs
%     dataTrue(:, ir) = dataTrue(:, ir) / norm(dataTrue(:, ir), 2);
% end

%% Prepare noisy data
sigma = 0.1;
noise = sigma * randn(size(dataTrue));
noisyData = dataTrue + noise;
trainData = noisyData;

%% Parameters for dictionary learning using sparse K-SVD
gain = 1.15;
trainBlockSize = 16;                        % for each dimension
trainBlockNum = 4096;                       % number of training blocks in the training set
trainIter = 10;
sigSpThres = sigma * trainBlockSize * gain; % pre-defined l2-norm error for BPDN
atomSpThres = 10;                           % a self-determind value to control the sparsity of matrix A

%% Base dictionary setting
% wavelet
nlevels_wavelet = [0, 0];	% Decomposition level, all 0 means wavelet
pfilter_wavelet = '9/7' ;	% Pyramidal filter
dfilter_wavelet = '9/7' ;	% Directional filter

%% Dictionary learning using sparse K-SVD
fprintf('------------------------------------------------------------\n');
fprintf('Dictionary Learning\n');

[vecTrainBlockCoeff, stb_wavelet] = pdfb2vec(pdfbdec(zeros(trainBlockSize, trainBlockSize), pfilter_wavelet, dfilter_wavelet, nlevels_wavelet));
initDict = speye(length(vecTrainBlockCoeff), length(vecTrainBlockCoeff));
baseSynOp = @(x) pdfb(x, stb_wavelet, pfilter_wavelet, dfilter_wavelet, nlevels_wavelet, trainBlockSize, trainBlockSize, 1);
baseAnaOp = @(x) pdfb(x, stb_wavelet, pfilter_wavelet, dfilter_wavelet, nlevels_wavelet, trainBlockSize, trainBlockSize, 2);
[learnedDict, Coeffs, err] = sparseKsvd(trainData, baseSynOp, baseAnaOp, ...
    initDict, trainIter, trainBlockSize, trainBlockNum, atomSpThres, sigSpThres, 'bpdn');

%% Denoising
fprintf('------------------------------------------------------------\n');
fprintf('Denoising\n');

cleanData = zeros(size(dataTrue));
% nSamples, nRecs
totalBlockNum = (nSamples - trainBlockSize + 1) * (nRecs - trainBlockSize + 1);
processedBlocks = 0;

for ibatch = 1:nRecs-trainBlockSize+1
    fprintf('Batch %d... ', ibatch);
    % the current batch of blocks
    blocks = im2colstep(noisyData(:, ibatch:ibatch+trainBlockSize-1), trainBlockSize * [1, 1], [1, 1]);
    cleanBlocks = zeros(size(blocks));
    
    blockCoeff = zeros(length(vecTrainBlockCoeff), nSamples - trainBlockSize + 1);
    for iblk = 1:nSamples - trainBlockSize + 1
        opts = spgSetParms('verbosity', 0, 'optTol', 1e-6);
        blockCoeff(:, iblk) = spg_bpdn(@(x, mode) learnedOp(x, baseSynOp, baseAnaOp, learnedDict, mode), blocks(:, iblk), sigSpThres, opts);
        % blockCoeff(:, iblk) = OMP({@(x) baseSynOp(learnedDict*x), @(x) learnedDict'*baseAnaOp(x)}, blocks(:, iblk), sigSpThres);
        cleanBlocks(:, iblk) =  learnedOp(blockCoeff(:, iblk), baseSynOp, baseAnaOp, learnedDict, 1);
    end
    
    cleanBatch = col2imstep(cleanBlocks, [nSamples, trainBlockSize], trainBlockSize * [1, 1], [1, 1]);
    cleanData(:,ibatch:ibatch+trainBlockSize-1) = cleanData(:,ibatch:ibatch+trainBlockSize-1) + cleanBatch;
    
    processedBlocks = processedBlocks + (nSamples - trainBlockSize + 1);
    fprintf('Processed %d blocks\n', processedBlocks);
end

% average the denoised and noisy signals
cnt = countcover(size(noisyData), trainBlockSize * [1, 1], [1, 1]);
cleanData = cleanData./cnt;

%% Plot figures and PSNR output
figure; imshow(dataTrue/max(dataTrue(:)));
title('Original Seismic Data');

figure; imshow(noisyData/max(noisyData(:)));
psnrNoisyData = 20*log10(max(dataTrue(:)) * sqrt(numel(noisyData)) / norm(dataTrue(:) - noisyData(:), 2));
title(sprintf('Noisy Seismic Data, PSNR = %.2fdB', psnrNoisyData));
fprintf('------------------------------------------------------------\n');
fprintf('Noisy Seismic Data, PSNR = %.2fdB\n', psnrNoisyData);

figure; imshow(cleanData/max(cleanData(:)));
psnrCleanData = 20*log10(max(dataTrue(:)) * sqrt(numel(cleanData)) / norm(dataTrue(:) - cleanData(:), 2));
title(sprintf('Denoised Seismic Data, PSNR = %.2fdB', psnrCleanData));
fprintf('------------------------------------------------------------\n');
fprintf('Denoised Seismic Data, PSNR = %.2fdB\n', psnrCleanData);