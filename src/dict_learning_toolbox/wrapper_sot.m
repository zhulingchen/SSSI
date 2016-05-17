function y = wrapper_sot(x, m, n, D, blkSize, blkSize_tlCorner, dirClass, dirRange, mode)
% WRAPPER_SOT is a wrapper function for sparse orthonormal transform (SOT)
% / inverse sparse orthonormal transform (iSOT)
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology


if (mode == 1) % iSOT
    x = reshape(x, prod(blkSize), []);
    
    % number of rows and columns interior blocks
    nIntBlockRows = floor((m - blkSize_tlCorner(1)) / blkSize(1));
    nIntBlockCols = floor((n - blkSize_tlCorner(2)) / blkSize(2));
    
    % size of the bottom-right corner block
    blkSize_brCorner = ([m, n] - blkSize_tlCorner) - [nIntBlockRows, nIntBlockCols] .* blkSize;
    
    nBlockRows = nIntBlockRows + 1 + (blkSize_brCorner(1) > 0);
    nBlockCols = nIntBlockCols + 1 + (blkSize_brCorner(2) > 0);
    
    x = mat2cell(x, prod(blkSize), ones(1, nBlockRows * nBlockCols));
    x = reshape(x, nBlockRows, nBlockCols);
    y = inverseSot(x, m, n, D, blkSize, blkSize_tlCorner, dirClass);
    y = y(:);
elseif (mode == 2) % SOT
    x = reshape(x, m, n);
    y = forwardSot(x, D, blkSize, blkSize_tlCorner, dirRange);
    y = [y{:}];
    y = y(:);
else
    error('Wrong mode!');
end