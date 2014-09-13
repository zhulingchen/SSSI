function y = sefilter2(x, f1, f2, extmod, shift)
% SEFILTER2   2D seperable filtering with extension handling
%
%       y = sefilter2(x, f1, f2, [extmod], [shift])
%
% Input:
%   x:      input image
%   f1, f2: 1-D filters in each dimension that make up a 2D seperable filter
%   extmod: [optional] extension mode (default is 'per')
%   shift:  [optional] specify the window over which the 
%       	convolution occurs. By default shift = [0; 0].
%
% Output:
%   y:      filtered image of the same size as the input image:
%           Y(z1,z2) = X(z1,z2)*F1(z1)*F2(z2)*z1^shift(1)*z2^shift(2)
%
% Note:
%   The origin of the filter f is assumed to be floor(size(f)/2) + 1.
%   Amount of shift should be no more than floor((size(f)-1)/2).
%   The output image has the same size with the input image.
%
% See also: EXTEND2, EFILTER2

if ~exist('extmod', 'var')
    extmod = 'per';
end

if ~exist('shift', 'var')
    shift = [0; 0];
end

% Make sure filter in a row vector
f1 = f1(:)';
f2 = f2(:)';

% Periodized extension
lf1 = (length(f1) - 1) / 2;
lf2 = (length(f2) - 1) / 2;

y = extend2(x, floor(lf1) + shift(1), ceil(lf1) - shift(1), ...
    floor(lf2) + shift(2), ceil(lf2) - shift(2), extmod);

% Seperable filter
y = conv2(f1, f2, y, 'valid');