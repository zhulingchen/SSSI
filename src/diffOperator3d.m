function diffData = diffOperator3d(data, c, d, dim)
% DIFFOPERATOR Performs higher-order approximation of staggered-grid finite difference
%
% data(n1, n2, n3)      input data
% n1                number of depth samples
% n2                number of horizontal samples
% c                 differentiation coefficients
% d                 distance per sample
% dim               dimension
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Entao Liu (liuentao@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

if (dim > 3)
    error('Dimension dim cannot be larger than 3!');
end

order = size(c, 1);
k = 2 * order - 1;
[n1, n2, n3] = size(data);

if (dim == 1)
    diffData = zeros(n1-k, n2, n3);
elseif dim == 2
    diffData = zeros(n1, n2-k, n3);
else
    diffData = zeros(n1, n2, n3-k);
end

for ii = 1:order
    if (dim == 1)
        idx1 = ((order+1):(n1-k+order)) + (ii-1);
        idx2 = ((order+1):(n1-k+order)) - ii;
        diffData = diffData + c(ii) * (data(idx1, :, :) - data(idx2, :, :)) / d;
    elseif (dim == 2)
        idx1 = ((order+1):(n2-k+order)) + (ii-1);
        idx2 = ((order+1):(n2-k+order)) - ii;
        diffData = diffData + c(ii) * (data(:, idx1, :) - data(:, idx2, :)) / d;
    elseif (dim == 3)
        idx1 = ((order+1):(n3-k+order)) + (ii-1);
        idx2 = ((order+1):(n3-k+order)) - ii;
        diffData = diffData + c(ii) * (data(:, :, idx1) - data(:, :, idx2)) / d;
    end
end
