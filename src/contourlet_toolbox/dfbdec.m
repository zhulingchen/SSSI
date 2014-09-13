function y = dfbdec(x, fname, n)
% DFBDEC   Directional Filterbank Decomposition
%
%	y = dfbdec(x, fname, n)
%
% Input:
%   x:      input image
%   fname:  filter name to be called by DFILTERS
%   n:      number of decomposition tree levels
%
% Output:
%   y:	    subband images in a cell vector of length 2^n
%
% Note:
%   This is the general version that works with any FIR filters
%      
% See also: DFBREC, FBDEC, DFILTERS

if (n ~= round(n)) | (n < 0)
    error('Number of decomposition levels must be a non-negative integer');
end

if n == 0
    % No decomposition, simply copy input to output
    y{1} = x;    
    return;
end

% Get the diamond-shaped filters
[h0, h1] = dfilters(fname, 'd');

% Fan filters for the first two levels
%   k0: filters the first dimension (row)
%   k1: filters the second dimension (column)
k0 = modulate2(h0, 'c');
k1 = modulate2(h1, 'c');    

% Tree-structured filter banks
if n == 1
    % Simplest case, one level
    [y{1}, y{2}] = fbdec(x, k0, k1, 'q', '1r', 'per');
    
else
    % For the cases that n >= 2
    % First level
    [x0, x1] = fbdec(x, k0, k1, 'q', '1r', 'per');
    
    % Second level
    y = cell(1, 4);
    [y{1}, y{2}] = fbdec(x0, k0, k1, 'q', '2c', 'qper_col');
    [y{3}, y{4}] = fbdec(x1, k0, k1, 'q', '2c', 'qper_col');
    
    % Fan filters from diamond filters
    [f0, f1] = ffilters(h0, h1);

    % Now expand the rest of the tree
    for l = 3:n
        % Allocate space for the new subband outputs
        y_old = y;    
        y = cell(1, 2^l);
	
        % The first half channels use R1 and R2
        for k = 1:2^(l-2)
            i = mod(k-1, 2) + 1;
            [y{2*k-1}, y{2*k}] = fbdec(y_old{k}, f0{i}, f1{i}, 'pq', i, 'per');
        end	
	
        % The second half channels use R3 and R4	
        for k = 2^(l-2)+1:2^(l-1)
            i = mod(k-1, 2) + 3;
            [y{2*k-1}, y{2*k}] = fbdec(y_old{k}, f0{i}, f1{i}, 'pq', i, 'per');
        end
    end
end

% Back sampling (so that the overal sampling is separable) 
% to enhance visualization
y = backsamp(y);

% Flip the order of the second half channels
y(2^(n-1)+1:end) = fliplr(y(2^(n-1)+1:end));