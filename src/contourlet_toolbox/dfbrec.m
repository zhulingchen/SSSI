function x = dfbrec(y, fname)
% DFBREC   Directional Filterbank Reconstruction
%
%   x = dfbrec(y, fname)
%
% Input:
%   y:	    subband images in a cell vector of length 2^n
%   fname:  filter name to be called by DFILTERS
%
% Output:
%   x:	    reconstructed image
%
% See also: DFBDEC, FBREC, DFILTERS

n = log2(length(y));

if (n ~= round(n)) | (n < 0)
    error('Number of reconstruction levels must be a non-negative integer');
end

if n == 0
    % Simply copy input to output
    x = y{1};
    return;
end

% Get the diamond-shaped filters
[h0, h1] = dfilters(fname, 'r');

% Fan filters for the first two levels
%   k0: filters the first dimension (row)
%   k1: filters the second dimension (column)
k0 = modulate2(h0, 'c');
k1 = modulate2(h1, 'c');

% Flip back the order of the second half channels
y(2^(n-1)+1:end) = fliplr(y(2^(n-1)+1:end));

% Undo backsampling
y = rebacksamp(y);

% Tree-structured filter banks
if n == 1
    % Simplest case, one level
    x = fbrec(y{1}, y{2}, k0, k1, 'q', '1r', 'per');
    
else
    % For the cases that n >= 2

    % Fan filters from diamond filters
    [f0, f1] = ffilters(h0, h1);

    % Recombine subband outputs to the next level
    for l = n:-1:3
        y_old = y;    
        y = cell(1, 2^(l-1));
	
        % The first half channels use R1 and R2
        for k = 1:2^(l-2)
            i = mod(k-1, 2) + 1;
            y{k} = fbrec(y_old{2*k-1}, y_old{2*k}, f0{i}, f1{i}, ...
                'pq', i, 'per');
        end
	
        % The second half channels use R3 and R4
        for k = 2^(l-2)+1:2^(l-1)
            i = mod(k-1, 2) + 3;
            y{k} = fbrec(y_old{2*k-1}, y_old{2*k}, f0{i}, f1{i}, ...
                'pq', i, 'per');
        end	
    end

    % Second level
    x0 = fbrec(y{1}, y{2}, k0, k1, 'q', '2c', 'qper_col');
    x1 = fbrec(y{3}, y{4}, k0, k1, 'q', '2c', 'qper_col');
    
    % First level
    x = fbrec(x0, x1, k0, k1, 'q', '1r', 'per');
end