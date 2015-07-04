%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Yue M. Lu and Minh N. Do
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	PyrNDRec_mm.m
%	
%   First Created: 10-11-05
%	Last Revision: 07-13-09
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rec = PyrNDRec_mm(subs, InD, Pyr_mode, smooth_func)

%   N-dimensional multiscale pyramid reconstruction - with multiple modes
%
%   INPUT:
%
%     subs: an L+1 by 1 cell array storing subbands from the coarsest scale
%           to the coarsest scale. subs{1} contains the lowpass subband.
%  
%     InD: "S", then the output subbands are given in the spatial domain.
%          "F", then the output subbands are given in the frequency domain.
%          This option can be used to get rid of the unnecessary (and
%          time-consuming) "ifftn-fftn" operations in the middle steps.
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
%   rec: the reconstructed signal in the spatial domain
%
%   See also:
%
%   PyrNDRec_mm.m


InD = upper(InD);

N = ndims(subs{1});

L = length(subs) - 1;

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

Lp = subs{1};
subs{1} = [];
if InD == 'S'
    Lp = fftn(Lp);
end
Lp = ccsym(Lp, N, 'c');

for n = 1 : L
    
    Hp = subs{n+1};
    subs{n+1} = [];
    
    if InD == 'S'
        Hp = fftn(Hp);
    end
    Hp = ccsym(Hp, N, 'c');
    Lp = PrySDrec_onestep(Lp, Hp, w_array(n), tbw_array(n), D_array(n), smooth_func);
     
end

clear Hp

Lp = ccsym(Lp, N, 'e');
rec = real(ifftn(Lp));

