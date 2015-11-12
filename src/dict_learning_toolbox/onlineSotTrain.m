function [Dout, X, Aout] = onlineSotTrain(Y, Din, t, Ain, nBlocks, option)
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


[szPatch, nAtoms] = size(Din);

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

D = Din;
mu = (t * nBlocks + 1 - nBlocks)/(t * nBlocks + 1);
for iter = 1:iterations
    % fix dictionary, hard thresholding coefficients
    X = D' * Y;
    X(abs(X) < lambda) = 0;
    
    cost_updateX = norm(Y - D * X, 'fro')^2 + nnz(X(:)) * lambda * lambda;
    
    % fix coefficients, update dictionary
    Z = mu * Ain + X * Y';
    [U, ~, V] = svd(Z);
    D = V * [U, zeros(nAtoms, szPatch - nAtoms)]';
    
    cost_updateD = norm(Y - D * X, 'fro')^2 + nnz(X(:)) * lambda * lambda;
    
    fprintf('iter = %d, cost_updateX = %f, cost_updateD = %f\n', iter, cost_updateX, cost_updateD);
    
    if ( abs(cost_updateD - cost_updateX) < tol )
        break;
    end
end

Dout = D;
Aout = mu * Ain + X * Y';