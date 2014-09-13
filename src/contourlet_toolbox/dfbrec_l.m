function x = dfbrec_l(y, f)
% DFBREC_L   Directional Filterbank Reconstruction using Ladder Structure
%
%	x = dfbrec_l(y, fname)
%
% Input:
%	y:	subband images in a cell vector of length 2^n
%	f:	filter in the ladder network structure,
%		can be a string naming a standard filter (see LDFILTER)
%
% Output:
%	x:	reconstructed image
%
% See also:	DFBDEC, FBREC, DFILTERS

n = log2(length(y));

if (n ~= round(n)) | (n < 0)
    error('Number of reconstruction levels must be a non-negative integer');
end

if n == 0
    % Simply copy input to output
    x = y{1};
    return;
end

% Ladder filter
if isstr(f)
    f = ldfilter(f);
end

% Flip back the order of the second half channels
y(2^(n-1)+1:end) = fliplr(y(2^(n-1)+1:end));

% Undo backsampling
y = rebacksamp(y);

% Tree-structured filter banks
if n == 1
    % Simplest case, one level
    x = fbrec_l(y{1}, y{2}, f, 'q', '1r', 'qper_col');
    
else
    % For the cases that n >= 2
    
    % Recombine subband outputs to the next level
    for l = n:-1:3
        y_old = y;
        y = cell(1, 2^(l-1));
        
        % The first half channels use R1 and R2
        for k = 1:2^(l-2)
            i = mod(k-1, 2) + 1;
            y{k} = fbrec_l(y_old{2*k}, y_old{2*k-1}, f, 'p', i, 'per');
        end

        % The second half channels use R3 and R4
        for k = 2^(l-2)+1:2^(l-1)
            i = mod(k-1, 2) + 3;
            y{k} = fbrec_l(y_old{2*k}, y_old{2*k-1}, f, 'p', i, 'per');
        end
    end
	
    % Second level
    x0 = fbrec_l(y{2}, y{1}, f, 'q', '2c', 'per');
    x1 = fbrec_l(y{4}, y{3}, f, 'q', '2c', 'per');

    % First level
    x = fbrec_l(x0, x1, f, 'q', '1r', 'qper_col');
end