function [p0, p1] = ppdec(x, type)
% PPDEC   Parallelogram Polyphase Decomposition
%
% 	[p0, p1] = ppdec(x, type)
%
% Input:
%	x:	input image
%	type:	one of {1, 2, 3, 4} for selecting sampling matrices:
%			P1 = [2, 0; 1, 1]
%			P2 = [2, 0; -1, 1]
%			P3 = [1, 1; 0, 2]
%			P4 = [1, -1; 0, 2]
%
% Output:
%	p0, p1:	two parallelogram polyphase components of the image
%
% Note:
%	These sampling matrices appear in the directional filterbank:
%		P1 = R1 * Q1
%		P2 = R2 * Q2
%		P3 = R3 * Q2
%		P4 = R4 * Q1
%	where R's are resampling matrices and Q's are quincunx matrices
%
% See also:	QPDEC

% Parallelogram polyphase decomposition by simplifying sampling matrices
% using the Smith decomposition of the quincunx matrices

switch type
    case 1	% P1 = R1 * Q1 = D1 * R3	
	p0 = resamp(x(1:2:end, :), 3);
	
	% R1 * [0; 1] = [1; 1]
	p1 = resamp(x(2:2:end, [2:end, 1]), 3);
		
    case 2	% P2 = R2 * Q2 = D1 * R4
	p0 = resamp(x(1:2:end, :), 4);
	
	% R2 * [1; 0] = [1; 0]
	p1 = resamp(x(2:2:end, :), 4);
	
    case 3	% P3 = R3 * Q2 = D2 * R1
	p0 = resamp(x(:, 1:2:end), 1);
	
	% R3 * [1; 0] = [1; 1]
	p1 = resamp(x([2:end, 1], 2:2:end), 1);
	
    case 4	% P4 = R4 * Q1 = D2 * R2
	p0 = resamp(x(:, 1:2:end), 2);
	
	% R4 * [0; 1] = [0; 1]
	p1 = resamp(x(:, 2:2:end), 2);

    otherwise
	error('Invalid argument type');
end
