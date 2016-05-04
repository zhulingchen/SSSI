function X = ifdct_usfft(C, is_real)

% ifdct_usfft - Inverse curvelet transform
%
% Input
%     C         Curvelet coefficients
%     isreal    Type of transform
%                   0: complex
%                   1: real
%
% Output
%     X         Image
%

[m,n] = size(C{end}{1});

if nargin < 2, is_real = 0; end;

nbscales = length(C);

allcurvelets = 0;

if(is_real)
    C = fdct_usfft_r2c(C);
end

% call mex function
X = ifdct_usfft_mex(m,n,nbscales, nbangles_coarse, allcurvelets, C);



