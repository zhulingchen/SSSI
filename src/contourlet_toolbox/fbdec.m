function [y0, y1] = fbdec(x, h0, h1, type1, type2, extmod)
% FBDEC   Two-channel 2D Filterbank Decomposition
%
%	[y0, y1] = fbdec(x, h0, h1, type1, type2, [extmod])
%
% Input:
%	x:	input image
%	h0, h1:	two decomposition 2D filters
%	type1:	'q', 'p' or 'pq' for selecting quincunx or parallelogram
%		downsampling matrix
%	type2:	second parameter for selecting the filterbank type
%		If type1 == 'q' then type2 is one of {'1r', '1c', '2r', '2c'}
%		If type1 == 'p' then type2 is one of {1, 2, 3, 4}
%			Those are specified in QDOWN and PDOWN
%		If type1 == 'pq' then same as 'p' except that
%		the paralellogram matrix is replaced by a combination 
%		of a  resampling and a quincunx matrices
%	extmod:	[optional] extension mode (default is 'per')
%
% Output:
%	y0, y1:	two result subband images
%
% Note:		This is the general implementation of 2D two-channel
% 		filterbank
%
% See also:	FBDEC_SP

if ~exist('extmod', 'var')
    extmod = 'per';
end

% For parallegoram filterbank using quincunx downsampling, resampling is
% applied before filtering
if type1 == 'pq'
    x = resamp(x, type2);
end

% Stagger sampling if filter is odd-size (in both dimensions)
if all(mod(size(h1), 2))
    shift = [-1; 0];
    
    % Account for the resampling matrix in the parallegoram case
    if type1 == 'p'
	R = {[1, 1; 0, 1], [1, -1; 0, 1], [1, 0; 1, 1], [1, 0; -1, 1]};	
	shift = R{type2} * shift;
    end        
    
else
    shift = [0; 0];
end

% Extend, filter and keep the original size
y0 = efilter2(x, h0, extmod);
y1 = efilter2(x, h1, extmod, shift);

% Downsampling
switch type1
    case 'q'
    	% Quincunx downsampling
    	y0 = qdown(y0, type2);
    	y1 = qdown(y1, type2);
	
    case 'p'
    	% Parallelogram downsampling
    	y0 = pdown(y0, type2);
    	y1 = pdown(y1, type2);
	
    case 'pq'
    	% Quincux downsampling using the equipvalent type
    	pqtype = {'1r', '2r', '2c', '1c'};

    	y0 = qdown(y0, pqtype{type2});
    	y1 = qdown(y1, pqtype{type2});
	
    otherwise
    	error('Invalid input type1');
end