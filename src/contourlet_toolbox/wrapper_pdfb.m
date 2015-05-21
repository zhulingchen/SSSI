function y = wrapper_pdfb(x, s, pfilt, dfilt, nlevs, nz, nx, mode)
% WRAPPER_PDFB Can be used as a wrapper function handle to run Pyramidal
% Directional Filter Bank (PDFB), or called Contourlet / inverse Contourlet
% transforms
%
% it must have the signature to be used as a function handle in spgl1
% package for l1-optimization
%
% A = @(x, mode) wrapper_pdfb(x, s, is_real, nz, nx, mode);
%
% y = A(x,mode)   if mode == 1 then y = A x  (y is m-by-1);
%                 if mode == 2 then y = A'x  (y is n-by-1).
%
% input arguments
% c             Contourlet coefficients in the form of vector
% s             structure of the Contourlet output
% mode          transform mode -- 1: inverse transform; 2: transform
%
% output arguments
% X             inverse transformed image in the form of vector
%
% See also:	PDFBDEC, PDFBREC, VEC2PDFB, PDFB2VEC
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology


if (mode == 1) % inverse Contourlet transform
    x = vec2pdfb(x, s);
    y = pdfbrec(x, pfilt, dfilt);
    y = y(:);
elseif (mode == 2) % Contourlet transform
    x = reshape(x, nz, nx);
    y = pdfbdec(x, pfilt, dfilt, nlevs);
    y = pdfb2vec(y);
else
    error('Wrong mode!');
end