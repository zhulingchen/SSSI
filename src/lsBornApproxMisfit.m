function [value, grad] = lsBornApproxMisfit(dm, w, fs, dataDeltaFreq, greenFreqForShotSet, greenFreqForRecSet)
% LSBORNAPPROXMISFIT Calculates the least-squares misfit function with
% respect to the perturbation model dm based on the Born approximation
%
% value = 1/2 * (L(dm) - delta_d)' * (L(dm) - delta_d)
% grad = real(L'(L(dm) - delta_d))
% where L is the forward modelling operator based on the Born approximation
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology


nLength = length(dm);
nw = length(w);
nShots = size(dataDeltaFreq, 2);

% value of the cost function
value = 0;
% gradient of the cost function
grad = zeros(nLength, 1);

% update the velocity model with least-squares
parfor iw = 1:nw
    
    % fprintf('Processing f = %fHz ... ', w(iw)/(2*pi));
    % tic;
    
    greenFreqForShot = greenFreqForShotSet{iw};
    greenFreqForRec = greenFreqForRecSet{iw};
    
    fieldScatter = w(iw)^2 * fs(iw) * (repmat(dm, 1, nShots) .* greenFreqForShot).' * greenFreqForRec;
    bias = fieldScatter.' - dataDeltaFreq(:, :, iw);
    value = value + 1/2 * norm(bias, 'fro')^2;
    
    grad = grad + w(iw)^2 * fs(iw) * sum(greenFreqForShot .* (greenFreqForRec * conj(bias)), 2);
    
    % timePerFreq = toc;
    % fprintf('elapsed time = %fs\n', timePerFreq);
    
end

grad = real(grad);

% fprintf('error function value = %f\n', value);