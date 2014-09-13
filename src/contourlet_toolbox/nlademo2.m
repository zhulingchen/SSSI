function nlademo2
% NLADEMO2  Nonlinear approximation demo using only the finest scale

% Test image
im = imread('peppers.png');
im = double(im) / 256;          % image range = [0, 1]
im = smthborder(im, 32);        % smooth the borders to avoid border effect

% Parameters
pfilt = '9-7';
dfilt = '9-7';
nlevs = [5];

% Wavevelet coefficients at the finest scale
wc = pdfbdec(im, pfilt, dfilt, [0]);

% Contourlet coefficients at the finest scale
cc = pdfbdec(im, pfilt, dfilt, nlevs);

scale = 1;

wave_nla = cell(1, 4);
pdfb_nla = cell(1, 4);

% Nonlinear approximation at the finest level
for k = 1:4,
    wave_nla{k} = pdfbrec(pdfb_tr(wc, scale, 0, 4^k), pfilt, dfilt);
    pdfb_nla{k} = pdfbrec(pdfb_tr(cc, scale, 0, 4^k), pfilt, dfilt);
end

range = [0, 0.002];

figure(1)
clf;
set(gcf, 'Name', 'Wavelets')
for k = 1:4,
    subplot(2,2,k),
    imagesc(wave_nla{k}, range), axis image, axis off
    str = sprintf('M = %d', 4^k);
    title(str);
end
colormap('gray(256)')

figure(2)
clf;
set(gcf, 'Name', 'Contourlets')
for k = 1:4,
    subplot(2,2,k),
    imagesc(pdfb_nla{k}, range), axis image, axis off
    str = sprintf('M = %d', 4^k);
    title(str);
end
colormap('gray(256)')