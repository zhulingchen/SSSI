function diffData = diffOperator(data, c, d, dim)
% DIFFOPERATOR Performs higher-order approximation of staggered-grid finite difference
%
% data(n1, n2)      input data (2d) or
% data(n1, n2, n3)	input data (3d) or
% n1                number of depth samples
% n2                number of horizontal samples
% c                 differentiation coefficients
% d                 distance per sample
% dim               dimension
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% It is a Matlab interface of its background C/Mex function to replace its
% Matlab function.
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

diffData = diffOperator_mex(data, c, d, dim);