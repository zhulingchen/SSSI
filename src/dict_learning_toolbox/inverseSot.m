function x = inverseSot(y, m, n, D, blkSize, blkSize_tlCorner, dirClass)
% INVERSESOT inverse sparse orthonormal transform (iSOT)
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology


x = zeros(m, n);

if (~all((blkSize - blkSize_tlCorner) >= 0))
    error('Corner block size is larger than interior block size!');
end

% number of rows and columns interior blocks
nIntBlockRows = floor((m - blkSize_tlCorner(1)) / blkSize(1));
nIntBlockCols = floor((n - blkSize_tlCorner(2)) / blkSize(2));

% size of the bottom-right corner block
blkSize_brCorner = (size(x) - blkSize_tlCorner) - [nIntBlockRows, nIntBlockCols] .* blkSize;

% process top-left corner block
if (iscell(D))
    recBlock = D{dirClass(1, 1)} * y{1, 1};
else
    recBlock = D * y{1, 1};
end
recBlock = reshape(recBlock, blkSize);
x(1:blkSize_tlCorner(1), 1:blkSize_tlCorner(2)) = x(1:blkSize_tlCorner(1), 1:blkSize_tlCorner(2)) + recBlock(1:blkSize_tlCorner(1), 1:blkSize_tlCorner(2));

% process top-right corner block
if (blkSize_brCorner(2) > 0)
    if (iscell(D))
        recBlock = D{dirClass(1, end)} * y{1, end};
    else
        recBlock = D * y{1, end};
    end
    recBlock = reshape(recBlock, blkSize);
    x(1:blkSize_tlCorner(1), end-blkSize_brCorner(2)+1:end) = x(1:blkSize_tlCorner(1), end-blkSize_brCorner(2)+1:end) + recBlock(1:blkSize_tlCorner(1), end-blkSize_brCorner(2)+1:end);
end

% process bottom-left corner block
if (blkSize_brCorner(1) > 0)
    if (iscell(D))
        recBlock = D{dirClass(end, 1)} * y{end, 1};
    else
        recBlock = D * y{end, 1};
    end
    recBlock = reshape(recBlock, blkSize);
    x(end-blkSize_brCorner(1)+1:end, 1:blkSize_tlCorner(2)) = x(end-blkSize_brCorner(1)+1:end, 1:blkSize_tlCorner(2)) + recBlock(end-blkSize_brCorner(1)+1:end, 1:blkSize_tlCorner(2));
end

% process bottom-right corner block
if (blkSize_brCorner(1) > 0 && blkSize_brCorner(2) > 0)
    if (iscell(D))
        recBlock = D{dirClass(end, end)} * y{end, end};
    else
        recBlock = D * y{end, end};
    end
    recBlock = reshape(recBlock, blkSize);
    x(end-blkSize_brCorner(1)+1:end, end-blkSize_brCorner(2)+1:end) = x(end-blkSize_brCorner(1)+1:end, end-blkSize_brCorner(2)+1:end) + recBlock(end-blkSize_brCorner(1)+1:end, end-blkSize_brCorner(2)+1:end);
end

% process top boundary blocks
for idxBlockCol = 1:nIntBlockCols
    idxCol = ((idxBlockCol-1)*blkSize(2)+1:idxBlockCol*blkSize(2)) + blkSize_tlCorner(2);
    if (iscell(D))
        recBlock = D{dirClass(1, idxBlockCol + 1)} * y{1, idxBlockCol + 1};
    else
        recBlock = D * y{1, idxBlockCol + 1};
    end
    recBlock = reshape(recBlock, blkSize);
    x(1:blkSize_tlCorner(1), idxCol) = x(1:blkSize_tlCorner(1), idxCol) + recBlock(1:blkSize_tlCorner(1), :);
end

% process left boundary blocks
for idxBlockRow = 1:nIntBlockRows
    idxRow = ((idxBlockRow-1)*blkSize(1)+1:idxBlockRow*blkSize(1)) + blkSize_tlCorner(1);
    if (iscell(D))
        recBlock = D{dirClass(idxBlockRow + 1, 1)} * y{idxBlockRow + 1, 1};
    else
        recBlock = D * y{idxBlockRow + 1, 1};
    end
    recBlock = reshape(recBlock, blkSize);
    x(idxRow, 1:blkSize_tlCorner(2)) = x(idxRow, 1:blkSize_tlCorner(2)) + recBlock(:, 1:blkSize_tlCorner(2));
end

% process right boundary blocks
if (blkSize_brCorner(2) > 0)
    for idxBlockRow = 1:nIntBlockRows
        idxRow = ((idxBlockRow-1)*blkSize(1)+1:idxBlockRow*blkSize(1)) + blkSize_tlCorner(1);
        if (iscell(D))
            recBlock = D{dirClass(idxBlockRow + 1, end)} * y{idxBlockRow + 1, end};
        else
            recBlock = D * y{idxBlockRow + 1, end};
        end
        recBlock = reshape(recBlock, blkSize);
        x(idxRow, end-blkSize_brCorner(2)+1:end) = x(idxRow, end-blkSize_brCorner(2)+1:end) + recBlock(:, end-blkSize_brCorner(2)+1:end);
    end
end

% process bottom boundary blocks
if (blkSize_brCorner(1) > 0)
    for idxBlockCol = 1:nIntBlockCols
        idxCol = ((idxBlockCol-1)*blkSize(2)+1:idxBlockCol*blkSize(2)) + blkSize_tlCorner(2);
        if (iscell(D))
            recBlock = D{dirClass(end, idxBlockCol + 1)} * y{end, idxBlockCol + 1};
        else
            recBlock = D * y{end, idxBlockCol + 1};
        end
        recBlock = reshape(recBlock, blkSize);
        x(end-blkSize_brCorner(1)+1:end, idxCol) = x(end-blkSize_brCorner(1)+1:end, idxCol) + recBlock(end-blkSize_brCorner(1)+1:end, :);
    end
end

% process interior blocks
for idxBlockCol = 1:nIntBlockCols
    
    idxCol = ((idxBlockCol-1)*blkSize(2)+1:idxBlockCol*blkSize(2)) + blkSize_tlCorner(2);
    % process block by block in the current column
    for idxBlockRow = 1:nIntBlockRows
        idxRow = ((idxBlockRow-1)*blkSize(1)+1:idxBlockRow*blkSize(1)) + blkSize_tlCorner(1);
        % inverse transform
        if (iscell(D))
            recBlock = D{dirClass(idxBlockRow + 1, idxBlockCol + 1)} * y{idxBlockRow + 1, idxBlockCol + 1};
        else
            recBlock = D * y{idxBlockRow + 1, idxBlockCol + 1};
        end
        recBlock = reshape(recBlock, blkSize);
        % tile the block back
        x(idxRow, idxCol) = x(idxRow, idxCol) + recBlock;
    end
    
end