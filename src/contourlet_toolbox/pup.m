function y = pup(x, type, phase)
% PUP   Parallelogram Upsampling
%
% 	y = pup(x, type, [phase])
%
% Input:
%	x:	input image
%	type:	one of {1, 2, 3, 4} for selecting sampling matrices:
%			P1 = [2, 0; 1, 1]
%			P2 = [2, 0; -1, 1]
%			P3 = [1, 1; 0, 2]
%			P4 = [1, -1; 0, 2]
%	phase:	[optional] 0 or 1 to specify the phase of the input image as
%		zero- or one-polyphase	component, (default is 0)
%
% Output:
%	y:	parallelogram upsampled image
%
% Note:
%	These sampling matrices appear in the directional filterbank:
%		P1 = R1 * Q1
%		P2 = R2 * Q2
%		P3 = R3 * Q2
%		P4 = R4 * Q1
%	where R's are resampling matrices and Q's are quincunx matrices
%
% See also:	PPDEC

if ~exist('phase', 'var')
    phase = 0;
end

% Parallelogram polyphase decomposition by simplifying sampling matrices
% using the Smith decomposition of the quincunx matrices
%
% Note that R1 * R2 = R3 * R4 = I so for example,
% upsample by R1 is the same with down sample by R2.
% Also the order of upsampling operations is in the reserved order
% with the one of matrix multiplication.

[m, n] = size(x);

switch type
    case 1	% P1 = R1 * Q1 = D1 * R3	
	y = zeros(2*m, n);
	
	if phase == 0
	    y(1:2:end, :) = resamp(x, 4);
	else
	    y(2:2:end, [2:end, 1]) = resamp(x, 4);
	end
	
    case 2	% P2 = R2 * Q2 = D1 * R4
	y = zeros(2*m, n);
	
	if phase == 0
	    y(1:2:end, :) = resamp(x, 3);
	else
	    y(2:2:end, :) = resamp(x, 3);
	end
	
    case 3	% P3 = R3 * Q2 = D2 * R1
	y = zeros(m, 2*n);
	
	if phase == 0
	    y(:, 1:2:end) = resamp(x, 2);
	else
	    y([2:end, 1], 2:2:end) = resamp(x, 2);
	end	
	
    case 4	% P4 = R4 * Q1 = D2 * R2
	y = zeros(m, 2*n);
	
	if phase == 0
	    y(:, 1:2:end) = resamp(x, 1);
	else
	    y(:, 2:2:end) = resamp(x, 1);
	end	

    otherwise
	error('Invalid argument type');
end
