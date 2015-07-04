%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Yue M. Lu and Minh N. Do
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	ContourletSD_ne.m
%	
%   First Created: 06-09-08
%	Last Revision: 07-13-09
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Estimate the norm scaling factor for each contourlet coefficient

% size of the image: N-by-N 
N = 512;

smooth_func = @rcos;
nlev_SD = [2 2 3 4 5];
dfilt = 'pkva';

% Pyr_mode can be set to be 1, 1.5 or 2. Correspondingly, the redundancy 
% factors of the transform will be 2.33, 1.59, or 1.33, respectively.
Pyr_mode = 2;

if Pyr_mode == 1.5
    N = 480;
end

% The number of repeated experiments
niter = 200;

h = waitbar(0, 'Please wait ...');
for n = 1 : niter
    X = randn(N, N);
    
    Y = ContourletSDDec(X, nlev_SD, Pyr_mode, smooth_func, dfilt);
    
    if n == 1
        E = cell(size(Y));
        for s = 2:length(Y)
            E{s} = cell(size(Y{s}));
            for w=1:length(Y{s})
                E{s}{w} = Y{s}{w} .^ 2;
            end
        end
    else
        for s = 2 : length(Y)
            for w = 1 : length(Y{s})
                E{s}{w} = E{s}{w} + Y{s}{w} .^ 2;
            end
        end
    end
    
    waitbar(n / niter, h);
end
close(h);

for s = 2 : length(Y)
    for w = 1 : length(Y{s})
        E{s}{w} = sqrt(E{s}{w} / (niter - 1));
    end
end

% save the output in a mat file
pmode = Pyr_mode;
if pmode == 1.5
    pmode = 15;
end
eval(['save SDmm_' num2str(N) '_' num2str(pmode) ' E']);
