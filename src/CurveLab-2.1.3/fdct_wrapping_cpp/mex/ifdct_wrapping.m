function X = ifdct_wrapping(C, is_real, M, N)

% ifdct_wrapping - Inverse curvelet transform
%
% Input
%     C         Curvelet coefficients
%     isreal    Type of transform
%                   0: complex
%                   1: real
%
% Output
%     X         A double precision matrix
%
% See also ifdct_wrapping in the fdct_wrapping_matlab/ directory.

nbscales = length(C);
nbangles_coarse = length(C{2});
if length(C{end}) == 1, finest = 2; else finest = 1; end;

if nargin < 2, is_real = 0; end;
if nargin < 4,
    if finest == 1, error('Syntax: IFCT_wrapping(C,M,N) where the matrix to be recovered is M-by-N'); end;
end;

if(is_real)
    C = fdct_wrapping_r2c(C);
end

% call mex function
X = ifdct_wrapping_mex(M, N, nbscales, nbangles_coarse, finest, C);



