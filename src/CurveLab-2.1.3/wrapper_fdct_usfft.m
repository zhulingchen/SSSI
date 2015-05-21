function y = wrapper_fdct_usfft(x, s, is_real, nbscales, nz, nx, mode)
% WRAPPER_FDCT_USFFT can be used as a wrapper function handle to run
% Curvelet / inverse Curvelet transforms based on the USFFT version
%
% it must have the signature to be used as a function handle in spgl1
% package for l1-optimization
%
% A = @(x, mode) wrapper_fdct_usfft(x, s, is_real, nbscales, nbangles_coarse, nz, nx, mode);
%
% y = A(x,mode)   if mode == 1 then y = A x  (y is m-by-1);
%                 if mode == 2 then y = A'x  (y is n-by-1).
%
% input arguments
% c             Curvelet coefficients in the form of vector
% s             structure of the Curvelet output
% is_real       type of transform -- 0: complex; 1: real
% mode          transform mode -- 1: inverse transform; 2: transform
%
% output arguments
% X             inverse transformed image in the form of vector
%
% See also:	FDCT_USFFT, IFDCT_USFFT, VEC2CURVELET, CURVELET2VEC
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology


if (mode == 1) % inverse Curvelet transform
    x = vec2curvelet(x, s);
    if ~isunix
        y = real(ifdct_usfft(x, is_real));
    else
        y = real(ifdct_usfft(x, is_real, nbscales));
    end
    y = y(:);
elseif (mode == 2) % Curvelet transform
    x = reshape(x, nz, nx);
    y = fdct_usfft(x, is_real, nbscales);
    y = curvelet2vec(y);
else
    error('Wrong mode!');
end