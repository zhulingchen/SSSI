function C = fdct_usfft(X, is_real, nbscales, nbangles_coarse)

% fdct_usfft - Forward curvelet transform
%
% Input
%     X         Image
%     isreal    Type of transform
%                   0: complex
%                   1: real
% Output
%     C         Curvelet coefficients
%

[m,n] = size(X);

if nargin < 2, is_real = 0; end;
if nargin < 3, nbscales = floor(log2(min(m,n)))-3; end;
if nargin < 4, nbangles_coarse = 16; end;

allcurvelets = 0;

%call mex function
C = fdct_usfft_mex(m, n, nbscales, nbangles_coarse, allcurvelets, double(X));

if(is_real)
    C = fdct_usfft_c2r(C);
end
