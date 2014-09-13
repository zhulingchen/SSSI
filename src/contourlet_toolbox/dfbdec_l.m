function y = dfbdec_l(x, f, n)
% DFBDEC_L   Directional Filterbank Decomposition using Ladder Structure
%
%	y = dfbdec_l(x, f, n)
%
% Input:
%	x:	input image
%	f:	filter in the ladder network structure,
%		can be a string naming a standard filter (see LDFILTER)
%	n:	number of decomposition tree levels
%
% Output:
%	y:	subband images in a cell array (of size 2^n x 1)

if (n ~= round(n)) | (n < 0)
    error('Number of decomposition levels must be a non-negative integer');
end

if n == 0
    % No decomposition, simply copy input to output
    y{1} = x;    
    return;
end

% Ladder filter
if isstr(f)
    f = ldfilter(f);
end

% Tree-structured filter banks
if n == 1
    % Simplest case, one level
    [y{1}, y{2}] = fbdec_l(x, f, 'q', '1r', 'qper_col');

else
    % For the cases that n >= 2
    % First level
    [x0, x1] = fbdec_l(x, f, 'q', '1r', 'qper_col');

    % Second level
    y = cell(1, 4);
    [y{2}, y{1}] = fbdec_l(x0, f, 'q', '2c', 'per');
    [y{4}, y{3}] = fbdec_l(x1, f, 'q', '2c', 'per');

    % Now expand the rest of the tree
    for l = 3:n
        % Allocate space for the new subband outputs
        y_old = y;
        y = cell(1, 2^l);

        % The first half channels use R1 and R2
        for k = 1:2^(l-2)
            i = mod(k-1, 2) + 1;
            [y{2*k}, y{2*k-1}] = fbdec_l(y_old{k}, f, 'p', i, 'per');
        end

        % The second half channels use R3 and R4        
        for k = 2^(l-2)+1:2^(l-1)
            i = mod(k-1, 2) + 3;
            [y{2*k}, y{2*k-1}] = fbdec_l(y_old{k}, f, 'p', i, 'per');
        end
    end
end

% Backsampling
y = backsamp(y);

% Flip the order of the second half channels
y(2^(n-1)+1:end) = fliplr(y(2^(n-1)+1:end));