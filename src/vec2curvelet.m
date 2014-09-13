function y = vec2curvelet(c, s)
% VEC2CURVELET   Convert the vector form to the output structure of the
% curvelet
%
%       y = vec2curvelet(c, s)
%
% Input:
%   c:  1-D vector that contains all curvelet coefficients
%   s:  structure of curvelet output
%
% Output:
%   y:  curvelet coefficients in cell vector format that can be used in ifdct_wrapping
%
% See also:	CURVELET2VEC, IFDCT_WRAPPING
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology


% Copy the coefficients from c to y according to the structure s
n = s(end, 1);      % number of pyramidal layers
y = cell(1, n);

% Variable that keep the current position
pos = 0;

% Used for row index of s
ind = 0;

for l = 1:n
    % Number of directional subbands in this layer
    nd = length(find(s(:, 1) == l));

    y{l} = cell(1, nd);
    
    for d = 1:nd
        % Size of this subband
        p = s(ind + d, 3);
        q = s(ind + d, 4);
        ss = p * q;
        
        y{l}{d} = reshape(c(pos+[1:ss]), [p, q]);
        pos = pos + ss;
    end
    
    ind = ind + nd;
end