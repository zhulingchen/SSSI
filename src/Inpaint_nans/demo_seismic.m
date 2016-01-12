close all;
clear;
clc;

load('cmp_cube_m2752.mat');
[nSamples, nTraces] = size(dataTrue);

% Remove mean and normalization
minData = min(dataTrue, [], 1);
maxData = max(dataTrue, [], 1);
% meanData = mean(dataTrue, 1);
% dataTrue = bsxfun(@times, bsxfun(@minus, dataTrue, meanData), 1./abs(maxData - minData));
dataTrue = bsxfun(@times, dataTrue, 1./abs(maxData - minData));
dataTrue = bsxfun(@times, dataTrue, 1./max(abs(dataTrue), [], 1)); % scale to [-1, 1]

% Prepare noisy data
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

% Prepare missing data
nNullTraces = round(0.2 * nTraces);
idxNullTraces = sort(randperm(nTraces, nNullTraces));
idxFullTraces = setdiff(1:nTraces, idxNullTraces);
noisyData(:, idxNullTraces) = NaN;

% inpainting
interpData = inpaint_nans(noisyData, 4);

figure;
subplot(1, 3, 1); imagesc(dataTrue); title('Original'); colormap(gray);
subplot(1, 3, 2); imagesc(noisyData); title('Noisy (Missing)'); colormap(gray);
subplot(1, 3, 3); imagesc(interpData); title('Noisy (Inpainted)'); colormap(gray);