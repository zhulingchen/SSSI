function Vnew = extBoundary(V, nBoundary, dim)
%
% EXTBOUNDARY extends layer for the 2-d or 3-d model on the bottom, left wall,
% right wall, front wall, and back wall (not the top!!), for the purpose of apply ABC.
%
% input argument
% V                 velocity model matrix
% nBoundary         thickness of the absorbing boundary
% dim               dimension (2 or 3)
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Entao Liu (liuentao@gmail.com), Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

if (nargin == 2)
    dim = 2;
end

if (dim == 2)
    
    % extension of left and right boundary
    Vnew = [repmat(V(:, 1), 1, nBoundary), V, repmat(V(:, end), 1, nBoundary)];
    
    % extension of bottom boundary
    Vnew(end+1:end+nBoundary, :) = repmat(Vnew(end,:), nBoundary, 1);
    
elseif (dim == 3)
    
    % extension of left and right boundary
    Vnew = [repmat(V(:, 1, :), [1, nBoundary, 1]), V, repmat(V(:, end, :), [1, nBoundary, 1])];
    
    % extension of front and rear boundary
    Vnew = cat(3, repmat(Vnew(:, :, 1), [1, 1, nBoundary]), Vnew, repmat(Vnew(:, :, end), [1, 1, nBoundary]));
    
    % extension of bottom boundary
    Vnew(end+1:end+nBoundary, :, :) = repmat(Vnew(end, :, :), [nBoundary, 1, 1]);
    
else
    error('d can only be 2 or 3!');
end

