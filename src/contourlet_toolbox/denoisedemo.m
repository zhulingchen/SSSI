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
im = imread('lena.png');
im = double(im) / 256;

% Generate noisy image. 
sig = std(im(:));
sigma = sig / rho;
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
range = [0, 1];

subplot(2,2,1), imagesc(im(41:168, 181:308), range); axis image off
set(gca, 'FontSize', 8);
title('Original Image', 'FontSize', 10);

subplot(2,2,2), imagesc(nim(41:168, 181:308), range); axis image off
set(gca, 'FontSize', 8);
title(sprintf('Noisy Image (SNR = %.2f dB)', ...
              SNR(im, nim)), 'FontSize', 10);

subplot(2,2,3), imagesc(wim(41:168, 181:308), range); axis image off
set(gca, 'FontSize', 8);
title(sprintf('Denoise using Wavelets (SNR = %.2f dB)', ...
              SNR(im, wim)), 'FontSize', 10);

subplot(2,2,4), imagesc(cim(41:168, 181:308), range); axis image off
set(gca, 'FontSize', 8);
title(sprintf('Denoise using Contourlets (SNR = %.2f dB)', ...
              SNR(im, cim)), 'FontSize', 10);