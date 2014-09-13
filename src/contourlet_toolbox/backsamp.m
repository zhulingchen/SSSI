function y = backsamp(y)
% BACKSAMP    Backsampling the subband images of the directional filter bank
%
%   y = backsamp(y)
%
% Input and output are cell vector of dyadic length
%
% This function is called at the end of the DFBDEC to obtain subband images 
% with overall sampling as diagonal matrices
%
% See also: DFBDEC

% Number of decomposition tree levels
n = log2(length(y));

if (n ~= round(n)) | (n < 1)
    error('Input must be a cell vector of dyadic length');
end

if n == 1
    % One level, the decomposition filterbank shoud be Q1r
    % Undo the last resampling (Q1r = R2 * D1 * R3)
    for k = 1:2
	y{k} = resamp(y{k}, 4);
	y{k}(:, 1:2:end) = resamp(y{k}(:, 1:2:end), 1);
	y{k}(:, 2:2:end) = resamp(y{k}(:, 2:2:end), 1);	
    end    
    
elseif n > 2    
    N = 2^(n-1);
    
    for k = 1:2^(n-2)
	shift = 2*k - (2^(n-2) + 1);
	
	% The first half channels
	y{2*k-1} = resamp(y{2*k-1}, 3, shift);
	y{2*k} = resamp(y{2*k}, 3, shift);
	
	% The second half channels
	y{2*k-1+N} = resamp(y{2*k-1+N}, 1, shift);
	y{2*k+N} = resamp(y{2*k+N}, 1, shift);	
    end
end