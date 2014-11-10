function flag=iscomplex(matrix)
%
% flag=iscomplex(matrix)
%
% returns 1 if any elements of the input matrix are complex. 
% returns 0 otherwise. This is exactly one line:
% flag=any(any(imag(matrix)));
%
% G.F. Margrave October 1993
flag = any(any(imag(matrix)));
