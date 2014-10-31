function denoisedemo
% DENOISEDEMO   Denoise demo
% Compare the denoise performance of wavelet and contourlet transforms
% Note: Noise standard deviation estimation of PDFB (function pdfb_nest) 
% can take a while...

% Parameters
pfilt = '9-7';
dfilt = 'pkva';
nlevs = [0, 0, 4, 4, 5];    % Number of levels for DFB at each pyramidal level
th = 3;                     % lead to 3*sigma threshold denoising
rho = 3;                    % noise level

% Test image: the usual suspect...
im = imread('barbara.png');
im = double(im);

% Generate noisy image. 
sig = std(im(:));
sigma = 0.1 * 255;%sig / rho;
nim = im + sigma * randn(size(im));


%%%%% Wavelet denoising %%%%%
% Wavelet transform using PDFB with zero number of level for DFB
y = pdfbdec(nim, pfilt, dfilt, zeros(length(nlevs), 1));
[c, s] = pdfb2vec(y);

% Threshold (typically 3*sigma)
wth = th * sigma;
c = c .* (abs(c) > wth);

% Reconstruction
y = vec2pdfb(c, s);
wim = pdfbrec(y, pfilt, dfilt);


%%%%% Contourlet Denoising %%%%%
% Contourlet transform
y = pdfbdec(nim, pfilt, dfilt, nlevs);
[c, s] = pdfb2vec(y);
c = c.';

% Threshold
% Require to estimate the noise standard deviation in the PDFB domain first 
% since PDFB is not an orthogonal transform
nvar = pdfb_nest(size(im,1), size(im, 2), pfilt, dfilt, nlevs);

cth = th * sigma * sqrt(nvar);
% cth = (4/3) * th * sigma * sqrt(nvar);

% Slightly different thresholds for the finest scale
fs = s(end, 1);
fssize = sum(prod(s(find(s(:, 1) == fs), 3:4), 2));
cth(end-fssize+1:end) = (4/3) * cth(end-fssize+1:end);

c = c .* (abs(c) > cth);

% Reconstruction
y = vec2pdfb(c, s);
cim = pdfbrec(y, pfilt, dfilt);


%%%%% Plot: Only the hat!
range = [0, 255];

% (41:168, 181:308)

MSE_noisy = sum(sum((im - nim).^2))/numel(im);
PSNR_noisy = 20*log10(255/sqrt(MSE_noisy));
MSE_waveletRestored = sum(sum((im - wim).^2))/numel(im);
PSNR_waveletRestored = 20*log10(255/sqrt(MSE_waveletRestored));
MSE_contourletRestored = sum(sum((im - cim).^2))/numel(im);
PSNR_contourletRestored = 20*log10(255/sqrt(MSE_contourletRestored));

subplot(2,2,1), imagesc(im, range); colormap gray; axis('image');
set(gca, 'FontSize', 8);
title('Original Image', 'FontSize', 10);

subplot(2,2,2), imagesc(nim, range); colormap gray; axis('image');
set(gca, 'FontSize', 8);
title(sprintf('Noisy Image (SNR = %.2f dB)', ...
              PSNR_noisy), 'FontSize', 10);

subplot(2,2,3), imagesc(wim, range); colormap gray; axis('image');
set(gca, 'FontSize', 8);
title(sprintf('Denoise using Wavelets (SNR = %.2f dB)', ...
              PSNR_waveletRestored), 'FontSize', 10);

subplot(2,2,4), imagesc(cim, range); colormap gray; axis('image');
set(gca, 'FontSize', 8);
title(sprintf('Denoise using Contourlets (SNR = %.2f dB)', ...
              PSNR_contourletRestored), 'FontSize', 10);