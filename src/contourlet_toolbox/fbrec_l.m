function x = fbrec_l(y0, y1, f, type1, type2, extmod)
% FBREC_L   Two-channel 2D Filterbank Reconstruction using Ladder Structure 
%
%	x = fbrec_l(y0, y1, f, type1, type2, [extmod])
%
% Input:
%	y0, y1:	two input subband images
%	f:	filter in the ladder network structure
%	type1:	'q' or 'p' for selecting quincunx or parallelogram
%		downsampling matrix
%	type2:	second parameter for selecting the filterbank type
%		If type1 == 'q' then type2 is one of {'1r', '1c', '2r', '2c'}
%			({2, 3, 1, 4} can also be used as equivalent)
%		If type1 == 'p' then type2 is one of {1, 2, 3, 4}
%		Those are specified in QPDEC and PPDEC
%       extmod: [optional] extension mode (default is 'per')
%		This refers to polyphase components.
%
% Output:
%	x:	reconstructed image
%
% Note:		This is also called the lifting scheme	
%
% See also:	FBDEC_L

% Modulate f
f(1:2:end) = -f(1:2:end);

if ~exist('extmod', 'var')
    extmod = 'per';
end

% Ladder network structure
p1 = (-1 / sqrt(2)) * (y1 + sefilter2(y0, f, f, extmod));
p0 = sqrt(2) * y0 + sefilter2(p1, f, f, extmod, [1, 1]);

% Polyphase reconstruction
switch lower(type1(1))
    case 'q'
        % Quincunx polyphase reconstruction
        x = qprec(p0, p1, type2);
	
    case 'p'
        % Parallelogram polyphase reconstruction
        x = pprec(p0, p1, type2);
        
    otherwise
        error('Invalid argument type1');
end