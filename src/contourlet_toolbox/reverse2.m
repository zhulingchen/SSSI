function r = reverse2(x)
% REVERSE2   Reverse order of elements in 2-d signal

if ndims(x)~=2, error('X must be a 2-D matrix.'); end

r = x(end:-1:1, end:-1:1);
