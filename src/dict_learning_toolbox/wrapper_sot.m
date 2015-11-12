function y = wrapper_sot(x, D, blkSize, nz, nx, dirClass, dirRange, mode)
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
    nBlockRows = floor(nz / blkSize(1));
    nBlockCols = floor(nx / blkSize(2));
    x = mat2cell(x, prod(blkSize), ones(1, nBlockRows * nBlockCols));
    x = reshape(x, nBlockRows, nBlockCols);
    y = inverseSot(x, D, blkSize, dirClass);
    y = y(:);
elseif (mode == 2) % SOT
    x = reshape(x, nz, nx);
    y = forwardSot(x, D, blkSize, dirRange);
    y = [y{:}];
    y = y(:);
else
    error('Wrong mode!');
end