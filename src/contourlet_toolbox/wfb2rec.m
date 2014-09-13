function x = wfb2rec(x_LL, x_LH, x_HL, x_HH, h, g)
% WFB2REC   2-D Wavelet Filter Bank Decomposition
%
%       x = wfb2rec(x_LL, x_LH, x_HL, x_HH, h, g)
%
% Input:
%   x_LL, x_LH, x_HL, x_HH:   Four 2-D wavelet subbands
%   h, g:   lowpass analysis and synthesis wavelet filters
%
% Output:
%   x:      reconstructed image

% Make sure filter in a row vector
h = h(:)';
g = g(:)';

g0 = g;
len_g0 = length(g0);
ext_g0 = floor((len_g0 - 1) / 2);

% Highpass synthesis filter: G1(z) = -z H0(-z)
len_g1 = length(h);
c = floor((len_g1 + 1) / 2);
g1 = (-1) * h .* (-1) .^ ([1:len_g1] - c);
ext_g1 = len_g1 - (c + 1);

% Get the output image size
[height, width] = size(x_LL);
x_B = zeros(height * 2, width);
x_B(1:2:end, :) = x_LL;

% Column-wise filtering
x_L = rowfiltering(x_B', g0, ext_g0)';
x_B(1:2:end, :) = x_LH;
x_L = x_L + rowfiltering(x_B', g1, ext_g1)';

x_B(1:2:end, :) = x_HL;
x_H = rowfiltering(x_B', g0, ext_g0)';
x_B(1:2:end, :) = x_HH;
x_H = x_H + rowfiltering(x_B', g1, ext_g1)';

% Row-wise filtering
x_B = zeros(2*height, 2*width);
x_B(:, 1:2:end) = x_L;
x = rowfiltering(x_B, g0, ext_g0);
x_B(:, 1:2:end) = x_H;
x = x + rowfiltering(x_B, g1, ext_g1);


% Internal function: Row-wise filtering with border handling 
function y = rowfiltering(x, f, ext1)
ext2 = length(f) - ext1 - 1;
x = [x(:, end-ext1+1:end) x x(:, 1:ext2)];
y = conv2(x, f, 'valid');