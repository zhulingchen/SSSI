function [c, s] = curvelet2vec(y)
% CURVELET2VEC   Convert the output of the curvelet transform into a vector form
%
%       [c, s] = curvelet2vec(y)
%
% Input:
%   y:  an output of the curvelet transform
%
% Output:
%   c:  1-D vector that contains all curvelet transform coefficients
%   s:  structure of the curvelet output, which is a four-column matrix.  Each row
%       of s corresponds to one subband y{l}{d} from y, in which the first two
%       entries are layer index l and direction index d and the last two
%       entries record the size of y{l}{d}.
%
% See also:	FDCT_WRAPPING, VEC2CURVELET
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology


n = length(y);

% Used for row index of s
ind = 0;

for l = 1:n
    nd = length(y{l});
    
    for d = 1:nd
        s(ind + d, :) = [l, d, size(y{l}{d})];
    end
    
    ind = ind + nd;
end

% The total number of directional coefficients
nc = sum(prod(s(:, 3:4), 2));
if (mod(nc, 2))
    nc = nc + 1; % to avoid odd vector length and cause problem in optimization
end

% Assign the coefficients to the vector c
c = zeros(nc, 1);

% Variable that keep the current position
pos = 0;

% Bandpass subbands
for l = 1:n    
    for d = 1:length(y{l})
        ss = prod(size(y{l}{d}));
        c(pos+[1:ss]) = y{l}{d}(:);
        pos = pos + ss;
    end
end