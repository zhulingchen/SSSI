function y = forwardSot(x, D)
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
blkSize = sqrt(sqrt(numel(D)));
nBlockRows = floor(m / blkSize);
nBlockCols = floor(n / blkSize);
y = cell(nBlockRows, nBlockCols);

for idxBlockCol = 1:nBlockCols
    
    idxCol = (idxBlockCol-1)*blkSize+1:idxBlockCol*blkSize;
    % process block by block in the current column
    for idxBlockRow = 1:nBlockRows
        idxRow = (idxBlockRow-1)*blkSize+1:idxBlockRow*blkSize;
        block = x(idxRow, idxCol);
        blockCoeff = D' * block(:);
        y{idxBlockRow, idxBlockCol} = blockCoeff(:);
    end
    
end