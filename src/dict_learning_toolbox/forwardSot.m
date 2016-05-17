function [y, dirClass] = forwardSot(x, D, blkSize, blkSize_tlCorner, dirRange)
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

if (~all((blkSize - blkSize_tlCorner) >= 0))
    error('Corner block size is larger than interior block size!');
end

% number of rows and columns interior blocks
nIntBlockRows = floor((m - blkSize_tlCorner(1)) / blkSize(1));
nIntBlockCols = floor((n - blkSize_tlCorner(2)) / blkSize(2));

% size of the bottom-right corner block
blkSize_brCorner = (size(x) - blkSize_tlCorner) - [nIntBlockRows, nIntBlockCols] .* blkSize;

y = cell(nIntBlockRows + 1 + (blkSize_brCorner(1) > 0), nIntBlockCols + 1 + (blkSize_brCorner(2) > 0));
dirClass = zeros(nIntBlockRows + 1 + (blkSize_brCorner(1) > 0), nIntBlockCols + 1 + (blkSize_brCorner(2) > 0));

% process top-left corner block
block = zeros(blkSize);
block(1:blkSize_tlCorner(1), 1:blkSize_tlCorner(2)) = x(1:blkSize_tlCorner(1), 1:blkSize_tlCorner(2));
if (iscell(D))
    % classify the block
    dirClass(1, 1) = getPatchGradClass(block, dirRange);
    % forward transform
    blockCoeff = D{dirClass(1, 1)}' * block(:);
else
    blockCoeff = D' * block(:);
end
y{1, 1} = blockCoeff(:);

% process top-right corner block
if (blkSize_brCorner(2) > 0)
    block = zeros(blkSize);
    block(1:blkSize_tlCorner(1), end-blkSize_brCorner(2)+1:end) = x(1:blkSize_tlCorner(1), end-blkSize_brCorner(2)+1:end);
    if (iscell(D))
        % classify the block
        dirClass(1, end) = getPatchGradClass(block, dirRange);
        % forward transform
        blockCoeff = D{dirClass(1, end)}' * block(:);
    else
        blockCoeff = D' * block(:);
    end
    y{1, end} = blockCoeff(:);
end

% process bottom-left corner block
if (blkSize_brCorner(1) > 0)
    block = zeros(blkSize);
    block(end-blkSize_brCorner(1)+1:end, 1:blkSize_tlCorner(2)) = x(end-blkSize_brCorner(1)+1:end, 1:blkSize_tlCorner(2));
    if (iscell(D))
        % classify the block
        dirClass(end, 1) = getPatchGradClass(block, dirRange);
        % forward transform
        blockCoeff = D{dirClass(end, 1)}' * block(:);
    else
        blockCoeff = D' * block(:);
    end
    y{end, 1} = blockCoeff(:);
end

% process bottom-right corner block
if (blkSize_brCorner(1) > 0 && blkSize_brCorner(2) > 0)
    block = zeros(blkSize);
    block(end-blkSize_brCorner(1)+1:end, end-blkSize_brCorner(2)+1:end) = x(end-blkSize_brCorner(1)+1:end, end-blkSize_brCorner(2)+1:end);
    if (iscell(D))
        % classify the block
        dirClass(end, end) = getPatchGradClass(block, dirRange);
        % forward transform
        blockCoeff = D{dirClass(end, end)}' * block(:);
    else
        blockCoeff = D' * block(:);
    end
    y{end, end} = blockCoeff(:);
end

% process top boundary blocks
for idxBlockCol = 1:nIntBlockCols
    idxCol = ((idxBlockCol-1)*blkSize(2)+1:idxBlockCol*blkSize(2)) + blkSize_tlCorner(2);
    block = zeros(blkSize);
    block(1:blkSize_tlCorner(1), :) = x(1:blkSize_tlCorner(1), idxCol);
    if (iscell(D))
        % classify the block
        dirClass(1, idxBlockCol + 1) = getPatchGradClass(block, dirRange);
        % forward transform
        blockCoeff = D{dirClass(1, idxBlockCol + 1)}' * block(:);
    else
        blockCoeff = D' * block(:);
    end
    y{1, idxBlockCol + 1} = blockCoeff(:);
end

% process left boundary blocks
for idxBlockRow = 1:nIntBlockRows
    idxRow = ((idxBlockRow-1)*blkSize(1)+1:idxBlockRow*blkSize(1)) + blkSize_tlCorner(1);
    block = zeros(blkSize);
    block(:, 1:blkSize_tlCorner(2)) = x(idxRow, 1:blkSize_tlCorner(2));
    if (iscell(D))
        % classify the block
        dirClass(idxBlockRow + 1, 1) = getPatchGradClass(block, dirRange);
        % forward transform
        blockCoeff = D{dirClass(idxBlockRow + 1, 1)}' * block(:);
    else
        blockCoeff = D' * block(:);
    end
    y{idxBlockRow + 1, 1} = blockCoeff(:);
end

% process right boundary blocks
if (blkSize_brCorner(2) > 0)
    for idxBlockRow = 1:nIntBlockRows
        idxRow = ((idxBlockRow-1)*blkSize(1)+1:idxBlockRow*blkSize(1)) + blkSize_tlCorner(1);
        block = zeros(blkSize);
        block(:, end-blkSize_brCorner(2)+1:end) = x(idxRow, end-blkSize_brCorner(2)+1:end);
        if (iscell(D))
            % classify the block
            dirClass(idxBlockRow + 1, end) = getPatchGradClass(block, dirRange);
            % forward transform
            blockCoeff = D{dirClass(idxBlockRow + 1, end)}' * block(:);
        else
            blockCoeff = D' * block(:);
        end
        y{idxBlockRow + 1, end} = blockCoeff(:);
    end
end

% process bottom boundary blocks
if (blkSize_brCorner(1) > 0)
    for idxBlockCol = 1:nIntBlockCols
        idxCol = ((idxBlockCol-1)*blkSize(2)+1:idxBlockCol*blkSize(2)) + blkSize_tlCorner(2);
        block = zeros(blkSize);
        block(end-blkSize_brCorner(1)+1:end, :) = x(end-blkSize_brCorner(1)+1:end, idxCol);
        if (iscell(D))
            % classify the block
            dirClass(end, idxBlockCol + 1) = getPatchGradClass(block, dirRange);
            % forward transform
            blockCoeff = D{dirClass(end, idxBlockCol + 1)}' * block(:);
        else
            blockCoeff = D' * block(:);
        end
        y{end, idxBlockCol + 1} = blockCoeff(:);
    end
end

% process interior blocks
for idxBlockCol = 1:nIntBlockCols
    
    idxCol = ((idxBlockCol-1)*blkSize(2)+1:idxBlockCol*blkSize(2)) + blkSize_tlCorner(2);
    % process block by block in the current column
    for idxBlockRow = 1:nIntBlockRows
        idxRow = ((idxBlockRow-1)*blkSize(1)+1:idxBlockRow*blkSize(1)) + blkSize_tlCorner(1);
        % extract the block
        block = x(idxRow, idxCol);
        if (iscell(D))
            % classify the block
            dirClass(idxBlockRow + 1, idxBlockCol + 1) = getPatchGradClass(block, dirRange);
            % forward transform
            blockCoeff = D{dirClass(idxBlockRow + 1, idxBlockCol + 1)}' * block(:);
        else
            blockCoeff = D' * block(:);
        end
        y{idxBlockRow + 1, idxBlockCol + 1} = blockCoeff(:);
    end
    
end