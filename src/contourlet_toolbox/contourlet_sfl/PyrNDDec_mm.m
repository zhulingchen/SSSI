%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Yue M. Lu and Minh N. Do
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	PyrNDDec_mm.m
%	
%   First Created: 10-11-05
%	Last Revision: 07-13-09
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subs = PyrNDDec_mm(X, OutD, L, Pyr_mode, smooth_func)

%   N-dimensional multiscale pyramid decomposition - with multiple modes
%
%   INPUT:
%
%     X: input signal in the spatial domain
%  
%     OutD: "S", then the output subbands are given in the spatial domain.
%           "F", then the output subbands are given in the frequency domain.
%           This option can be used to get rid of the unnecessary (and
%           time-consuming) "ifftn-fftn" operations in the middle steps.
%     L: level of decomposition
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
%   OUTPUT:
%
%   subs: an L+1 by 1 cell array storing subbands from the coarsest (i.e.
%   lowpass) to the finest scales.
%
%   See also:
%
%   PyrNDRec_mm.m

OutD = upper(OutD);

% The dimensionality of the input signal
N = ndims(X);

switch Pyr_mode
    case {1}
        % the cutoff frequencies at each scale
        w_array = [0.25 * ones(1, L - 1) 0.5];
        % the transition bandwidths at each scale
        tbw_array = [1/12 * ones(1, L - 1) 1/6];
        % the downsampling factor at each scale
        % no downsampling at the finest scale        
        D_array = [2 * ones(1, L - 1) 1];
        
    case {1.5}
        % the cutoff frequencies at each scale        
        w_array = [3/8 * ones(1, L - 1) 0.5];
        % the transition bandwidths at each scale        
        tbw_array = [1/9 * ones(1, L - 1) 1/7];
        % the downsampling factor at each scale
        % the lowpass channel at the first level of decomposition is
        % downsampled by a factor of 1.5 along each dimension.
        D_array = [2 * ones(1, L - 1) 1.5];
        
    case {2}
        % the cutoff frequencies at each scale                
        w_array = 1 / 3 * ones(1, L);
        % the transition bandwidths at each scale
        tbw_array = 1 / 7 * ones(1, L);
        % the downsampling factor at each scale        
        D_array = 2 * ones(1, L);
        
    otherwise
        error('Unsupported Pyr mode.');
end

X = fftn(X);
% We assume real-valued input signal X. Half of its Fourier
% coefficients can be removed due to conjugate symmetry.
X = ccsym(X, N, 'c');

subs = cell(L+1, 1);

for n = L : -1: 1
    % One level of the pyramid decomposition
    [Lp, Hp] = PrySDdec_onestep(X, w_array(n), tbw_array(n), D_array(n), smooth_func);
    
    X = Lp;
    
    Hp = ccsym(Hp, N, 'e'); 
    if OutD == 'S'
        % Go back to the spatial domain
        subs{n+1} = real(ifftn(Hp));
    else
        % Retain the frequency domain results
        subs{n+1} = Hp;
    end
    
    clear Lp Hp;
    
end

X = ccsym(X, N, 'e');
if OutD == 'S'
    subs{1} = real(ifftn(X));
else
    subs{1} = X;
end

        
        