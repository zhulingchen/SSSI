function C = fdct_wrapping(X, is_real, finest, nbscales, nbangles_coarse)

% fdct_wrapping - Forward curvelet transform
%
% Inputs
%     X         a double precision matrix
%     isreal    Type of transform
%                   0: complex
%                   1: real
%
% Output
%     C         Curvelet coefficients
%
% See also fdct_wrapping.m in the fdct_wrapping_matlab/ directory.

[M, N] = size(X);

if nargin < 2, is_real = 0; end;
if nargin < 3, nbscales = floor(log2(min(M, N)))-3; end;
if nargin < 4, nbangles_coarse = 16; end;

%call mex function
C = fdct_wrapping_mex(M, N, nbscales, nbangles_coarse, finest, double(X));

if(is_real)
    C = fdct_wrapping_c2r(C);
end
