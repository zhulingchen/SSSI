function y = smthborder(x, n)
% SMTHBORDER  Smooth the borders of a signal or image
%       y = smthborder(x, n)
%
% Input:
%   x:      the input signal or image
%   n:      number of samples near the border that will be smoothed
%
% Output:
%   y:      output image
%
% Note: This function provides a simple way to avoid border effect.

% Hamming window of size 2N
w = 0.54 - 0.46*cos(2*pi*(0:2*n-1)/(2*n-1));

if ndims(x) == 1
    W = ones(size(x));
    W(1:n) = w(1:n);
    W(end-n+1:end) = w(end-n+1:end);
    
    y = W .* x;
    
elseif ndims(x) == 2
    [n1, n2] = size(x);
    
    W1 = ones(n1, 1);
    W1(1:n) = w(1:n);
    W1(end-n+1:end) = w(n+1:end);
    
    y = W1(:, ones(1, n2)) .* x;
    
    W2 = ones(1, n2);
    W2(1:n) = w(1:n)';
    W2(end-n+1:end) = w(n+1:end)';
    
    y = W2(ones(n1, 1), :) .* y;
    
else
    error('First input must be a signal or image');
end