function diffData = diffOperator(data, c, d, dim)
% DIFFOPERATOR Performs higher-order approximation of staggered-grid finite difference
%
% data(n1, n2)      input data (2d) or
% data(n1, n2, n3)	input data (3d) or
% n1                number of depth samples
% n2                number of horizontal samples
% c                 differentiation coefficients
% d                 distance per sample
% dim               dimension
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

if nargin < 4
    dim = 1;
end

Ndims = ndims(data);

order = size(c, 1);
k = 2 * order - 1;
if (Ndims <= 2) % 2D case
    [n1, n2] = size(data);
    if (dim == 1)
        diffData = zeros(n1-k, n2);
        for ii = 1:order
            idx1 = ((order+1):(n1-k+order)) + (ii-1);
            idx2 = ((order+1):(n1-k+order)) - ii;
            diffData = diffData + c(ii) * (data(idx1, :) - data(idx2, :)) / d;
        end
    else % dim == 2
        diffData = zeros(n1, n2-k);
        for ii = 1:order
            idx1 = ((order+1):(n2-k+order)) + (ii-1);
            idx2 = ((order+1):(n2-k+order)) - ii;
            diffData = diffData + c(ii) * (data(:, idx1) - data(:, idx2)) / d;
        end
    end
else           % 3D case
    [n1, n2, n3] = size(data);
    if (dim == 1)
        diffData = zeros(n1-k, n2, n3);
        for ii = 1:order
            idx1 = ((order+1):(n1-k+order)) + (ii-1);
            idx2 = ((order+1):(n1-k+order)) - ii;
            diffData = diffData + c(ii) * (data(idx1, :, :) - data(idx2, :, :)) / d;
        end
    elseif (dim == 2)
        diffData = zeros(n1, n2-k, n3);
        for ii = 1:order
            idx1 = ((order+1):(n2-k+order)) + (ii-1);
            idx2 = ((order+1):(n2-k+order)) - ii;
            diffData = diffData + c(ii) * (data(:, idx1, :) - data(:, idx2, :)) / d;
        end
    else % dim == 3
        diffData = zeros(n1, n2, n3-k);
        for ii = 1:order
            idx1 = ((order+1):(n3-k+order)) + (ii-1);
            idx2 = ((order+1):(n3-k+order)) - ii;
            diffData = diffData + c(ii) * (data(:, :, idx1) - data(:, :, idx2)) / d;
        end
    end
end

