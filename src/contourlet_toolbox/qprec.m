function x = qprec(p0, p1, type)
% QPREC   Quincunx Polyphase Reconstruction
%
% 	x = qprec(p0, p1, [type])
%
% Input:
%	p0, p1:	two qunincunx polyphase components of the image
%	type:	[optional] one of {'1r', '1c', '2r', '2c'}, default is '1r'
%		'1' and '2' for selecting the quincunx matrices:
%			Q1 = [1, -1; 1, 1] or Q2 = [1, 1; -1, 1]
%		'r' and 'c' for suppresing row or column		
%
% Output:
%	x:	reconstructed image
%
% Note:
%	Note that R1 * R2 = R3 * R4 = I so for example,
%	upsample by R1 is the same with down sample by R2	
% 
% See also:	QPDEC

if ~exist('type', 'var')
    type = '1r';
end

% Quincunx downsampling using the Smith decomposition:
%
%       Q1 = R2 * D1 * R3
%          = R3 * D2 * R2
% and,
%       Q2 = R1 * D1 * R4
%          = R4 * D2 * R1
%
% where D1 = [2, 0; 0, 1] and D2 = [1, 0; 0, 2].
% See RESAMP for the definition of the resampling matrices R's

[m, n] = size(p0);

switch type
    case {'1r'}		% Q1 = R2 * D1 * R3
	y = zeros(2*m, n);
	
	y(1:2:end, :) = resamp(p0, 4);
	y(2:2:end, [2:end, 1]) = resamp(p1, 4);
	
	x = resamp(y, 1);
    
    case {'1c'}		% Q1 = R3 * D2 * R2
	y = zeros(m, 2*n);
	
	y(:, 1:2:end) = resamp(p0, 1);
	y(:, 2:2:end) = resamp(p1, 1);
	
	x = resamp(y, 4);
    
    case {'2r'}		% Q2 = R1 * D1 * R4
	y = zeros(2*m, n);
	
	y(1:2:end, :) = resamp(p0, 3);
	y(2:2:end, :) = resamp(p1, 3);
	
	x = resamp(y, 2);
		
    case {'2c'}		% Q2 = R4 * D2 * R1
	y = zeros(m, 2*n);
	
	y(:, 1:2:end) = resamp(p0, 2);
	y([2:end, 1], 2:2:end) = resamp(p1, 2);
	
	x = resamp(y, 3);
	
    otherwise
	error('Invalid argument type');
end
