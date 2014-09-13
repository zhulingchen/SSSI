function x = pprec(p0, p1, type)
% PPREC   Parallelogram Polyphase Reconstruction
%
% 	x = pprec(p0, p1, type)
%
% Input:
%	p0, p1:	two parallelogram polyphase components of the image
%	type:	one of {1, 2, 3, 4} for selecting sampling matrices:
%			P1 = [2, 0; 1, 1]
%			P2 = [2, 0; -1, 1]
%			P3 = [1, 1; 0, 2]
%			P4 = [1, -1; 0, 2]
%
% Output:
%	x:	reconstructed image
%
% Note:
%	These sampling matrices appear in the directional filterbank:
%		P1 = R1 * Q1
%		P2 = R2 * Q2
%		P3 = R3 * Q2
%		P4 = R4 * Q1
%	where R's are resampling matrices and Q's are quincunx matrices
%
%	Also note that R1 * R2 = R3 * R4 = I so for example,
%	upsample by R1 is the same with down sample by R2	
%
% See also:	PPDEC

% Parallelogram polyphase decomposition by simplifying sampling matrices
% using the Smith decomposition of the quincunx matrices

[m, n] = size(p0);

switch type
    case 1	% P1 = R1 * Q1 = D1 * R3
	x = zeros(2*m, n);
	
	x(1:2:end, :) = resamp(p0, 4);
	x(2:2:end, [2:end, 1]) = resamp(p1, 4);
		
    case 2	% P2 = R2 * Q2 = D1 * R4
	x = zeros(2*m, n);
	
	x(1:2:end, :) = resamp(p0, 3);
	x(2:2:end, :) = resamp(p1, 3);
	
    case 3	% P3 = R3 * Q2 = D2 * R1
	x = zeros(m, 2*n);
	
	x(:, 1:2:end) = resamp(p0, 2);
	x([2:end, 1], 2:2:end) = resamp(p1, 2);
	
    case 4	% P4 = R4 * Q1 = D2 * R2
	x = zeros(m, 2*n);
	
	x(:, 1:2:end) = resamp(p0, 1);
	x(:, 2:2:end) = resamp(p1, 1);

    otherwise
	error('Invalid argument type');
end
