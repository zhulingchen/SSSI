%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Yue M. Lu and Minh N. Do
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	ContourletSDRec.m
%	
%   First Created: 10-13-05
%	Last Revision: 07-13-09
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = ContourletSDRec(y, Pyr_mode, smooth_func, dfilt)

%   ContourletSD Reconstruction
%
%   INPUT
%
%     y: a cell vector of length length(nlevs) + 1, where except y{1} is 
%        the lowpass subband, each cell corresponds to one pyramidal
%        level and is a cell vector that contains bandpass directional
%        subbands from the DFB at that level.
%
%     Pyr_mode: Decomposition modes, including the following options.
%               1:   do not downsample the lowpass subband at the first level
%                    of decomposition
%               1.5: downsampling the first lowpass subbabnd by a factor of
%                    1.5 along each dimension
%               2:   use a lowpass filter with 1/3 pi cutoff frequency at the
%                    first level of decomposition.
%
%     smooth_func: function handle to generate the filter for the pyramid
%                  decomposition
%
%     dfilt: filter name for the directional decomposition step
%
%   OUTPUT
%
%     x: the reconstructed image
%
%   See also:
%
%   ContourletSDDec.m


L = length(y) - 1;

for k = 2 : L+1
    
    % Reconstruct the bandpass image from DFB

    % Decide the method based on the filter name
    switch dfilt
        case {'pkva6', 'pkva8', 'pkva12', 'pkva'}
            % Use the ladder structure (much more efficient)
            y{k} = dfbrec_l(y{k}, dfilt);

        otherwise
            % General case
            y{k} = dfbrec(y{k}, dfilt);
    end
end

x = PyrNDRec_mm(y, 'S', Pyr_mode, smooth_func);    
