function ytr = pdfb_tr(y, s, d, ncoef)
% PDFB_TR   Retain the most significant coefficients at certain subbands
%
%	ytr = pdfb_tr(y, s, d, [ncoef])
%
% Input:
%   y:      output from PDFB
%   s:      scale index (1 is the finest); 0 for ALL scales
%   d:      direction index; 0 for ALL directions
%   ncoef:  [optional] number of most significant coefficients from the
%           specified subbands; default is ALL coefficients.
%
% Output
%   ytr:    truncated PDFB output

n = length(y);
ytr = cell(1, n);

% Lowpass subband
if (s == n) | (s == 0)
    ytr{1} = y{1};
else
    ytr{1} = zeros(size(y{1}));
end

for l = 2:n
    for k = 1:length(y{l})
	    if (s == (n + 1 - l) | (s == 0)) & ((d == k) | (d == 0))
	        ytr{l}{k} = y{l}{k};	    
	    else
	        ytr{l}{k} = zeros(size(y{l}{k}));
	    end
    end
end

if exist('ncoef', 'var')  % Only keep ncoef most significant coefficients
    % Convert the output into the vector format
    [c, s] = pdfb2vec(ytr);

    % Sort the coefficient in the order of energy.
    csort = sort(abs(c));
    csort = fliplr(csort);
    
    thresh = csort(min(ncoef, length(csort)));
    cc = c .* (abs(c) >= thresh);
     
    ytr = vec2pdfb(cc, s);    
end
