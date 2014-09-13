function C = fdct_wrapping(X, isreal)

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
  
  [m,n] = size(X);
  nbscales = floor(log2(min(m,n)))-3;
  nbangles_coarse = 16;
  allcurvelets = 0;
  
  %call mex function
  C = fdct_wrapping_mex(m,n,nbscales, nbangles_coarse, allcurvelets, double(X));
  
  if(isreal)
    C = fdct_wrapping_c2r(C);
  end
