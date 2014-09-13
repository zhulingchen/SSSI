function x = pdfbrec(y, pfilt, dfilt)
% PDFBREC   Pyramid Directional Filterbank Reconstruction
%
%	x = pdfbrec(y, pfilt, dfilt)
%
% Input:
%   y:	    a cell vector of length n+1, one for each layer of 
%       	subband images from DFB, y{1} is the low band image
%   pfilt:  filter name for the pyramid
%   dfilt:  filter name for the directional filter bank
%
% Output:
%   x:      reconstructed image
%
% See also: PFILTERS, DFILTERS, PDFBDEC

n = length(y) - 1;
if n <= 0
    x = y{1};
    
else
    % Recursive call to reconstruct the low band
    xlo = pdfbrec(y(1:end-1), pfilt, dfilt);
    
    % Get the pyramidal filters from the filter name
    [h, g] = pfilters(pfilt);
    
    % Process the detail subbands
    if length(y{end}) ~= 3
        % Reconstruct the bandpass image from DFB
        
        % Decide the method based on the filter name
        switch dfilt        
            case {'pkva6', 'pkva8', 'pkva12', 'pkva'}	
                % Use the ladder structure (much more efficient)
                xhi = dfbrec_l(y{end}, dfilt);
                
            otherwise	
                % General case
                xhi = dfbrec(y{end}, dfilt); 
        end
        
        x = lprec(xlo, xhi, h, g);
   
    else    
        % Special case: length(y{end}) == 3
        % Perform one-level 2-D critically sampled wavelet filter bank
        x = wfb2rec(xlo, y{end}{1}, y{end}{2}, y{end}{3}, h, g);
    end
end