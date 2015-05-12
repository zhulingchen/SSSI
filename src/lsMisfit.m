function [value, grad] = lsMisfit(m, w, fs, dataTrueFreq, nz, nx, xs, zs, xr, zr, nDiffOrder, nBoundary, dz, dx)
% LSMISFIT Calculates the least-squares misfit function with respect to the
% model m defined by the differences at the receiver positions between the
% recorded seismic data and the modeled seismic data for each
% source-receiver pair of the seismic survey
%
% value = 1/2 * (F(m) - d_obs)' * (F(m) - d_obs)
% grad = real(F'(F(m) - d_obs))
% where F is the forward modelling operator
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology


nLength = numel(m);
nw = length(w);
nShots = length(xs);
nRecs = length(xr);

% velocity model (slowness)
m = reshape(m, nz + nBoundary, nx + 2*nBoundary);

% value of the cost function
value = 0;
% gradient of the cost function
grad = zeros(nLength, 1);

% update the velocity model with least-squares
parfor iw = 1:nw
    
    % received true data for all shots in frequency domain for current frequency
    sourceFreq = zeros(nLength, nShots);
    sourceFreq((xs-1)*(nz+nBoundary)+zs, :) = eye(nShots, nShots);
    [~, greenFreqForShot] = freqCpmlFor2dAw(m, sourceFreq, w(iw), nDiffOrder, nBoundary, dz, dx);
    
    % Green's function for every receiver
    sourceFreq = zeros(nLength, nRecs);
    sourceFreq((xr-1)*(nz+nBoundary)+zr, :) = eye(nRecs, nRecs);
    [~, greenFreqForRec] = freqCpmlFor2dAw(m, sourceFreq, w(iw), nDiffOrder, nBoundary, dz, dx);
    
    % get received data on the receivers
    dataCal = fs(iw) * greenFreqForShot((xr-1)*(nz+nBoundary)+zr, :);
    bias = dataCal - dataTrueFreq(:, :, iw); % dataTrueFreq(1:nRecs, 1:nShots, 1:nw)
    value = value + 1/2 * norm(bias, 'fro')^2;
    
    grad = grad + w(iw)^2 * fs(iw) * sum(greenFreqForShot .* (greenFreqForRec * conj(bias)), 2);
    
end

grad = real(grad);