function y = fdct(x, s, is_real, nz, nx, mode)
% FDCT Can be used as a function handle to run Curvelet / inverse
% Curvelet transforms
%
% it must have the signature to be used as a function handle in spgl1
% package for l1-optimization
%
% A = @(x, mode) fdct_handle(x, s, is_real, nz, nx, mode);
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
% See also:	FDCT_WRAPPING, IFDCT_WRAPPING, VEC2CURVELET, CURVELET2VEC
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
    y = ifdct_wrapping(x, is_real);
    y = y(:);
elseif (mode == 2) % Curvelet transform
    x = reshape(x, nz, nx);
    y = fdct_wrapping(x, is_real);
    y = curvelet2vec(y);
else
    error('Wrong mode!');
end