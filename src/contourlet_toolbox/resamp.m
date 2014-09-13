function y = resamp(x, type, shift, extmod)
% RESAMP   Resampling in 2D filterbank
%
%	y = resamp(x, type, [shift, extmod])
%
% Input:
%	x:	input image
%	type:	one of {1, 2, 3, 4} (see note)
%	shift:	[optional] amount of shift (default is 1)
%       extmod: [optional] extension mode (default is 'per').
%		Other options are:
%
% Output:
%	y:	resampled image.
%
% Note:
%	The resampling matrices are:
%		R1 = [1, 1;  0, 1];
%		R2 = [1, -1; 0, 1];
%		R3 = [1, 0;  1, 1];
%		R4 = [1, 0; -1, 1];
%
%	For type 1 and type 2, the input image is extended (for example
%	periodically) along the vertical direction;
%	while for type 3 and type 4 the image is extended along the 
%	horizontal direction.
%
%	Calling resamp(x, type, n) which n is positive integer is equivalent
%	to repeatly calling resamp(x, type) n times.
%
%	Input shift can be negative so that resamp(x, 1, -1) is the same
%	with resamp(x, 2, 1)

if ~exist('shift', 'var')
    shift = 1;
end

if ~exist('extmod', 'var')
    extmod = 'per';
end

switch type
    case {1, 2}
	y = resampc(x, type, shift, extmod);
	
    case {3, 4}
	y = resampc(x.', type - 2, shift, extmod).';
	
    otherwise
	error('The second input (type) must be one of {1, 2, 3, 4}');	
end