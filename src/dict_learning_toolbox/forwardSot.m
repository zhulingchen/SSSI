function [y, dirClass] = forwardSot(x, D, blkSize, dirRange)
% FORWARDSOT forward sparse orthonormal transform (SOT)
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology


[m, n] = size(x);
nBlockRows = floor(m / blkSize(1));
nBlockCols = floor(n / blkSize(2));
y = cell(nBlockRows, nBlockCols);
dirClass = zeros(nBlockRows, nBlockCols);

for idxBlockCol = 1:nBlockCols
    
    idxCol = (idxBlockCol-1)*blkSize(2)+1:idxBlockCol*blkSize(2);
    % process block by block in the current column
    for idxBlockRow = 1:nBlockRows
        idxRow = (idxBlockRow-1)*blkSize(1)+1:idxBlockRow*blkSize(1);
        % extract the block
        block = x(idxRow, idxCol);
        if (iscell(D))
            % classify the block
            dirClass(idxBlockRow, idxBlockCol) = getPatchGradClass(block, dirRange);
            % forward transform
            blockCoeff = D{dirClass(idxBlockRow, idxBlockCol)}' * block(:);
        else
            blockCoeff = D' * block(:);
        end
        y{idxBlockRow, idxBlockCol} = blockCoeff(:);
    end
    
end