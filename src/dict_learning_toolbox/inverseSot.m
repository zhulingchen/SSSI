function x = inverseSot(y, D, blkSize, dirClass)
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
m = nBlockRows * blkSize(1);
n = nBlockCols * blkSize(2);
x = zeros(m, n);

for idxBlockCol = 1:nBlockCols
    
    idxCol = (idxBlockCol-1)*blkSize(2)+1:idxBlockCol*blkSize(2);
    % process block by block in the current column
    for idxBlockRow = 1:nBlockRows
        idxRow = (idxBlockRow-1)*blkSize(1)+1:idxBlockRow*blkSize(1);
        % inverse transform
        if (iscell(D))
            recBlock = D{dirClass(idxBlockRow, idxBlockCol)} * y{idxBlockRow, idxBlockCol};
        else
            recBlock = D * y{idxBlockRow, idxBlockCol};
        end
        recBlock = reshape(recBlock, blkSize);
        % tile the block back
        x(idxRow, idxCol) = x(idxRow, idxCol) + recBlock;
    end
    
end