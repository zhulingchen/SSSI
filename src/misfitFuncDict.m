function [value, grad] = misfitFuncDict(dcoeff, invDictMat, dictMat, w, fs, nShots, dataDeltaFreq, greenFreqForShotSet, greenFreqForRecSet)
% MISFITFUNCDICT Calculate the least-squares misfit function with respect
% to the coefficients of perturbation model dm under some transform with
% explicit dictionary
%
% value = 1/2 * (L(invDictMat * dcoeff) - delta_d)' * (L(invDictMat * dcoeff) - delta_d)
% grad = dictMat * (L'(L(invDictMat * dcoeff) = delta_d))
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

nw = length(w);

% model after inverse transform
dm = invDictMat * dcoeff;
nLength = length(dm);

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
    
    fieldScatter = w(iw)^2 * (-fs(iw)) * (repmat(dm, 1, nShots) .* greenFreqForShot).' * greenFreqForRec;
    bias = fieldScatter - dataDeltaFreq(:, :, iw).';
    value = value + 1/2 * norm(bias, 'fro')^2;
    
    grad = grad + w(iw)^2 * (-fs(iw)) * sum(greenFreqForShot .* (greenFreqForRec * conj(bias.')), 2);
    
    % timePerFreq = toc;
    % fprintf('elapsed time = %fs\n', timePerFreq);
    
end

grad = real(dictMat * grad);

% fprintf('error function value = %f\n', value);

end