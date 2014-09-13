function diffData = diffOperator(data, c, d, dim)
% DIFFOPERATOR Performs higher-order approximation of staggered-grid finite difference
%
% data(n1, n2)      input data
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

if (dim > 2)
    error('Dimension dim cannot be larger than 2!');
end

order = size(c, 1);
k = 2 * order - 1;
[n1, n2] = size(data);

if (dim == 1)
    diffData = zeros(n1-k, n2);
else % dim == 2
    diffData = zeros(n1, n2-k);
end

for ii = 1:order
    if (dim == 1)
        idx1 = ((order+1):(n1-k+order)) + (ii-1);
        idx2 = ((order+1):(n1-k+order)) - ii;
        diffData = diffData + c(ii) * (data(idx1, :) - data(idx2, :)) / d;
    else % dim == 2
        idx1 = ((order+1):(n2-k+order)) + (ii-1);
        idx2 = ((order+1):(n2-k+order)) - ii;
        diffData = diffData + c(ii) * (data(:, idx1) - data(:, idx2)) / d;
    end
end
