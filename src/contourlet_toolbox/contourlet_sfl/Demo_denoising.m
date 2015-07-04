%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Yue M. Lu and Minh N. Do
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Demo_denoising.m
%	
%   First Created: 09-02-05
%	Last Revision: 07-13-09
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Denoising using the ContourletSD transform by simple hardthreholding

dfilt = 'pkva';
nlev_SD = [2 2 3 4 5];
smooth_func = @rcos;

% Pyr_mode can take the value of 1, 1.5 or 2. It specifies which one of
% the three variants of ContourletSD to use. The main differences are the
% redundancies of the transforms, which are 2.33, 1.59, or 1.33,
% respectively.

Pyr_mode = 1; 
% redundancy = 2. Set Pyr_mode = 2 for a less redundant version of the
% transform.



X = imread('lena.png');
X = double(X);

sigma = 30; % noise intensity

% Get the noisy image
Xn = X + sigma * randn(size(X));

% load the pre-computed norm scaling factor for each subband. For images
% of sizes other than 512 * 512, run ContourletSDmm_ne.m again.
eval(['load SDmm_' num2str(size(Xn,1)) '_' num2str(Pyr_mode)]);

% Take the ContourletSD transform
Y = ContourletSDDec(Xn, nlev_SD, Pyr_mode, smooth_func, dfilt);

% The redundancy of the transform can be verified as follows.
dstr1 = whos('Y');
dstr2 = whos('Xn');
dstr1.bytes / dstr2.bytes

% Apply hard thresholding on coefficients
Yth = Y;
for m = 2:length(Y)
  thresh = 3*sigma + sigma * (m == length(Y));
  for k = 1:length(Y{m})
    Yth{m}{k} = Y{m}{k}.* (abs(Y{m}{k}) > thresh*E{m}{k});
  end
end

% ContourletSD reconstruction
Xd = ContourletSDRec(Yth, Pyr_mode, smooth_func, dfilt);

figure
imshow(Xd / 255);
title(['Denoising using the contourletSD transform. PSNR (after denoising) = ' num2str(PSNR(X, Xd)) ' dB']);



% Denoising using the original contourlet transform

% Parameters
pfilt = '9-7';
dfilt = 'pkva';

nlevs = [0, 0, 4, 4, 5];    % Number of levels for DFB at each pyramidal level
% nlevs can be changed to, e.g., [2 2 3 4 5] (the same one used by ContourletSD above), 
% but the results are similar.
th = 3;                     % lead to 3*sigma threshold denoising

%%%%% Contourlet Denoising %%%%%
% Contourlet transform
y = pdfbdec(Xn, pfilt, dfilt, nlevs);
[c, s] = pdfb2vec(y);

% Threshold
% Require to estimate the noise standard deviation in the PDFB domain first
% since PDFB is not an orthogonal transform
% This step can take a while.
nvar = pdfb_nest(size(X,1), size(X, 2), pfilt, dfilt, nlevs);

cth = th * sigma * sqrt(nvar);
% cth = (4/3) * th * sigma * sqrt(nvar);

% Slightly different thresholds for the finest scale
fs = s(end, 1);
fssize = sum(prod(s(find(s(:, 1) == fs), 3:4), 2));
cth(end-fssize+1:end) = (4/3) * cth(end-fssize+1:end);

c = c .* (abs(c) > cth(:));

% Reconstruction
y = vec2pdfb(c, s);
Ycontour = pdfbrec(y, pfilt, dfilt);

figure;
imshow(Ycontour/255)
title(['Denoising using the original contourlet transform. PSNR (after denoising) = ' num2str(PSNR(X, Ycontour)) ' dB']);
