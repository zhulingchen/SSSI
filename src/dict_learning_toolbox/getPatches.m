function [Yout, dirClass] = getPatches(Yin, blkSize, nBlocks, dirRange)
% GETPATCHES Extract patches from input image
% GETPATCHES extracts patches from an input image and identifies the
% gradient classification of each patch
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

% create block training data
idx = cell(1, ndims(Yin));
if (isscalar(nBlocks))
    [idx{:}] = reggrid(size(Yin)-blkSize+1, nBlocks, 'eqnum');
else
    for idim = 1:ndims(Yin)
        idx{idim} = reggrid(size(Yin, idim)-blkSize(idim)+1, nBlocks(idim), 'eqnum');
    end
end
Yout = sampgrid(Yin, blkSize, idx{:});
nBlocks = size(Yout, 2);
dirClass = zeros(nBlocks, 1);
for iblock = 1:nBlocks
    % normalization
    if (norm(Yout(:, iblock), 2))
        Yout(:, iblock) = Yout(:, iblock) / norm(Yout(:, iblock), 2);
    end
    % gradient class estimation
    dirClass(iblock) = getPatchGradClass(reshape(Yout(:, iblock), blkSize), dirRange);
end