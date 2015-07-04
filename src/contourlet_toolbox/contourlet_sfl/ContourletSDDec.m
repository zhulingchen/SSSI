%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Yue M. Lu and Minh N. Do
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	ContourletSDDec.m
%	
%   First Created: 10-13-05
%	Last Revision: 07-13-09
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = ContourletSDDec(x, nlevs, Pyr_mode, smooth_func, dfilt)

%   ContourletSD Decomposition
%
%   INPUT
%
%     x: the input image.
%
%     nlevs: vector of numbers of directional filter bank decomposition
%            levels at each pyramidal level (from coarse to fine scale).          
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
%     y: a cell vector of length length(nlevs) + 1, where except y{1} is 
%        the lowpass subband, each cell corresponds to one pyramidal
%        level and is a cell vector that contains bandpass directional
%        subbands from the DFB at that level.
%
%   See also:
%
%   ContourletSDRec.m

L = length(nlevs);

y = PyrNDDec_mm(x, 'S', L, Pyr_mode, smooth_func);

for k = 2 : L+1
    % DFB on the bandpass image
    switch dfilt        % Decide the method based on the filter name
        case {'pkva6', 'pkva8', 'pkva12', 'pkva'}
            % Use the ladder structure (whihc is much more efficient)
            y{k} = dfbdec_l(y{k}, dfilt, nlevs(k-1));

        otherwise
            % General case
            y{k} = dfbdec(y{k}, dfilt, nlevs(k-1));
    end
end