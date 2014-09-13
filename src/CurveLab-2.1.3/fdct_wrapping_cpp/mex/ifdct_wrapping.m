function X = ifdct_wrapping(C, isreal)

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
  nbscales = floor(log2(min(m,n)))-3;
  nbangles_coarse = 16;
  allcurvelets = 0;
  
  if(isreal)
    C = fdct_wrapping_r2c(C);
  end
  
  % call mex function
  X = ifdct_wrapping_mex(m,n,nbscales, nbangles_coarse, allcurvelets, C);
  

  
