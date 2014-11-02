function X = ifdct_wrapping(C, is_real, nbscales, nbangles_coarse)

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

[m,n] = size(C{end}{1});

if nargin < 2, is_real = 0; end;
if nargin < 3, nbscales = floor(log2(min(m,n)))-3; end;
if nargin < 4, nbangles_coarse = 16; end;

allcurvelets = 0;

if(is_real)
    C = fdct_wrapping_r2c(C);
end

% call mex function
X = ifdct_wrapping_mex(m,n,nbscales, nbangles_coarse, allcurvelets, C);



