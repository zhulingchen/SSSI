function x = inverseSot(y, D)
% INVERSESOT inverse sparse orthonormal transform (iSOT)
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology


[nBlockRows, nBlockCols] = size(y);
blkSize = sqrt(sqrt(numel(D)));
m = nBlockRows * blkSize;
n = nBlockCols * blkSize;
x = zeros(m, n);

for idxBlockCol = 1:nBlockCols
    
    idxCol = (idxBlockCol-1)*blkSize+1:idxBlockCol*blkSize;
    % process block by block in the current column
    for idxBlockRow = 1:nBlockRows
        idxRow = (idxBlockRow-1)*blkSize+1:idxBlockRow*blkSize;
        recBlock = D * y{idxBlockRow, idxBlockCol};
        recBlock = reshape(recBlock, blkSize, blkSize);
        x(idxRow, idxCol) = x(idxRow, idxCol) + recBlock;
    end
    
end