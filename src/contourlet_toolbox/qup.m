function y = qup(x, type, phase)
% QUP   Quincunx Upsampling
%
% 	y = qup(x, [type], [phase])
%
% Input:
%	x:	input image
%	type:	[optional] one of {'1r', '1c', '2r', '2c'} (default is '1r')
%		'1' or '2' for selecting the quincunx matrices:
%			Q1 = [1, -1; 1, 1] or Q2 = [1, 1; -1, 1]
%		'r' or 'c' for extending row or column		
%	phase:	[optional] 0 or 1 to specify the phase of the input image as
%		zero- or one-polyphase component, (default is 0)
%
% Output:
%	y:	qunincunx upsampled image
%
% See also:	QDOWN

if ~exist('type', 'var')
    type = '1r';
end

if ~exist('phase', 'var')
    phase = 0;
end

% Quincunx downsampling using the Smith decomposition:
%
%	Q1 = R2 * [2, 0; 0, 1] * R3
%	   = R3 * [1, 0; 0, 2] * R2
% and,
%	Q2 = R1 * [2, 0; 0, 1] * R4
%	   = R4 * [1, 0; 0, 2] * R1
%
% See RESAMP for the definition of those resampling matrices
%
% Note that R1 * R2 = R3 * R4 = I so for example,
% upsample by R1 is the same with down sample by R2.
% Also the order of upsampling operations is in the reserved order
% with the one of matrix multiplication.

[m, n] = size(x);

switch type
    case {'1r'}
        z = zeros(2*m, n);
	
	if phase == 0
	    z(1:2:end, :) = resamp(x, 4);
	else
	    z(2:2:end, [2:end, 1]) = resamp(x, 4);
	end	
	
	y = resamp(z, 1);
	
    case {'1c'}
        z = zeros(m, 2*n);

	if phase == 0
	    z(:, 1:2:end) = resamp(x, 1);
	else
	    z(:, 2:2:end) = resamp(x, 1);
	end

	y = resamp(z, 4);
    
    case {'2r'}
        z = zeros(2*m, n);
	
	if phase == 0
	    z(1:2:end, :) = resamp(x, 3);
	else
	    z(2:2:end, :) = resamp(x, 3);
	end
	
	y = resamp(z, 2);
		
    case {'2c'}
        z = zeros(m, 2*n);
	
	if phase == 0
	    z(:, 1:2:end) = resamp(x, 2);
	else
	    z([2:end, 1], 2:2:end) = resamp(x, 2);
	end

	y = resamp(z, 3);
	
    otherwise
	error('Invalid argument type');
end
