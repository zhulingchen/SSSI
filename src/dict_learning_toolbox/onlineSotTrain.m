function [Dnew, X, Anew] = onlineSotTrain(Y, Dold, t, Aold, blkSize, nBlocks, option)
% ONLINESOTTRAIN Online sparse orthonormal dictionary learning
% ONLINESOTTRAIN runs dictionary learning in an online manner by using the
% previous dictionary as a warm start
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

dim = ndims(Y);
if ( dim < 2 || dim > 3 )
    error('Only 2-D and 3-D signals are supported!');
end

[szPatch, nAtoms] = size(Dold);

% create block training data
idx = cell(dim, 1);
[idx{:}] = reggrid(size(Y)-blkSize+1, nBlocks, 'eqdist');
Y = sampgrid(Y, blkSize, idx{:});
nBlocks = size(Y, 2);
for iblock = 1:nBlocks
    Y(:, iblock) = Y(:, iblock) / norm(Y(:, iblock), 2);
end

% parameter setting
if (isfield(option, 'lambda'))
    lambda = option.lambda;
else
    lambda = 0.1;
end

if (isfield(option, 'iterations'))
    iterations = option.iterations;
else
    iterations = 100;
end

if (isfield(option, 'tol'))
    tol = option.tol;
else
    tol = 1e-3;
end

D = Dold;
mu = (t * nBlocks + 1 - nBlocks)/(t * nBlocks + 1);
for iter = 1:iterations
    % fix dictionary, hard thresholding coefficients
    X = D' * Y;
    X(abs(X) < lambda) = 0;
    
    cost_updateX = norm(Y - D*X, 'fro')^2 + nnz(X(:)) * lambda * lambda;
    
    % fix coefficients, update dictionary
    Z = mu * Aold + X * Y';
    [U, ~, V] = svd(Z);
    D = V * [U, zeros(nAtoms, szPatch - nAtoms)]';
    
    cost_updateD = norm(Y - D*X, 'fro')^2 + nnz(X(:)) * lambda * lambda;
    
    fprintf('iter = %d, cost_updateX = %f, cost_updateD = %f\n', iter, cost_updateX, cost_updateD);
    
    if ( abs(cost_updateD - cost_updateX) < tol )
        break;
    end
end

Dnew = D;
Anew = mu * Aold + X * Y';