function y = dup(x, step, phase)
% DUP   Diagonal Upsampling
%
%	y = dup(x, step, [phase])
%
% Input:
%	x:	input image
%	step:	upsampling factors for each dimension which should be a
%		2-vector
%	phase:	[optional] to specify the phase of the input image which
%		should be less than step, (default is [0, 0])
%		If phase == 'minimum', a minimum size of upsampled image
%		is returned
%
% Output:
%	y:	diagonal upsampled image
%
% See also:	DDOWN

if ~exist('phase', 'var')
    phase = [0, 0];
end

sx = size(x);

if lower(phase(1)) == 'm'
    y = zeros((sx - 1) .* step + 1);    
    y(1:step(1):end, 1:step(2):end) = x;
    
else
    y = zeros(sx .* step);
    y(1+phase(1):step(1):end, 1+phase(2):step(2):end) = x;
end
