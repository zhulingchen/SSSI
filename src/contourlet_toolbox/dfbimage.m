function im = dfbimage(y, gap, gridI)
% DFBIMAGE    Produce an image from the result subbands of DFB
%
%	im = dfbimage(y, [gap, gridI])
%
% Input:
%	y:	output from DFBDEC
%	gap:	gap (in pixels) between subbands
%	gridI:	intensity of the grid that fills in the gap
%
% Output:
%	im:	an image with all DFB subbands
%
% The subband images are positioned as follows 
% (for the cases of 4 and 8 subbands):
%
%     0   1              0   1
%              and       2   3
%     2   3            4 5 6 7

% Gap between subbands
if ~exist('gap', 'var')
    gap = 0;    
end

l = length(y);

% Intensity of the grid (default is white)
if ~exist('gridI', 'var')
    gridI = 0;			
    for k = 1:l
	    m = max(abs(y{k}(:)));
	    % m = Inf;
	
	    if m > gridI
	        gridI = m;
	    end
    end
    
    % gridI = gridI * 1.1;		% add extra 10% of intensity
end

% Add grid seperation if required
if gap > 0
    for k = 1:l
	    y{k}(1:gap,:) = gridI;
	    y{k}(:,1:gap) = gridI;
    end	
end

% Simple case, only 2 subbands
if l == 2
    im = [y{1}; y{2}];
    return;
end

% Assume that the first subband has "horizontal" shape
[m, n] = size(y{1});

% The image
im = zeros(l*m/2, 2*n);

% First half of subband images ("horizontal" ones)
for k = 1:(l/4)
    im([1:m] + (k-1)*m, :) = [y{2*k-1}, y{2*k}];
end

% Second half of subband images ("vertical" ones)

% The size of each of those subband        
% It must be that: p = l*m/4  and n = l*q/4
[p, q] = size(y{l/2+1});
			
for k = 1:(l/2)
    im(p+1:end, [1:q] + (k-1)*q) = y{(l/2)+k};
end

% Finally, grid line in bottom and left
% if gap > 0
%     im(end-gap+1:end, :) = gridI;
%     im(:, end-gap+1:end) = gridI;
% end