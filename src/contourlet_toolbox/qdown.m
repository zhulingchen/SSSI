function y = qdown(x, type, extmod, phase)
% QDOWN   Quincunx Downsampling
%
% 	y = qdown(x, [type], [extmod], [phase])
%
% Input:
%	x:	input image
%	type:	[optional] one of {'1r', '1c', '2r', '2c'} (default is '1r')
%		'1' or '2' for selecting the quincunx matrices:
%			Q1 = [1, -1; 1, 1] or Q2 = [1, 1; -1, 1] 
%		'r' or 'c' for suppresing row or column		
%	phase:	[optional] 0 or 1 for keeping the zero- or one-polyphase
%		component, (default is 0)
%
% Output:
%	y:	qunincunx downsampled image
%
% See also:	QPDEC

if ~exist('type', 'var')
    type = '1r';
end

if ~exist('phase', 'var')
    phase = 0;
end

% Quincunx downsampling using the Smith decomposition:
%	Q1 = R2 * [2, 0; 0, 1] * R3
%	   = R3 * [1, 0; 0, 2] * R2
% and,
%	Q2 = R1 * [2, 0; 0, 1] * R4
%	   = R4 * [1, 0; 0, 2] * R1
%
% See RESAMP for the definition of those resampling matrices

switch type
    case {'1r'}
	z = resamp(x, 2);
	
	if phase == 0
	    y = resamp(z(1:2:end, :), 3);
	else
	    y = resamp(z(2:2:end, [2:end, 1]), 3);
	end	
	
    case {'1c'}
	z = resamp(x, 3);

	if phase == 0
	    y = resamp(z(:, 1:2:end), 2);
	else
	    y = resamp(z(:, 2:2:end), 2);
	end
	
    case {'2r'}
	z = resamp(x, 1);
	
	if phase == 0
	    y = resamp(z(1:2:end, :), 4);
	else
	    y = resamp(z(2:2:end, :), 4);
	end
		
    case {'2c'}
	z = resamp(x, 4);
	
	if phase == 0
	    y = resamp(z(:, 1:2:end), 1);
	else
	    y = resamp(z([2:end, 1], 2:2:end), 1);
	end

    otherwise
	error('Invalid argument type');
end
