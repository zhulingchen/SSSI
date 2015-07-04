%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Yue M. Lu and Minh N. Do
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Demo_basis_images.m
%	
%   First Created: 09-02-05
%	Last Revision: 07-13-09
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Display two basis functions of the ContourletSD in the spatial and
%  frequency domain.


% Set the number of directional subbands at different scales (from coarse
% to fine)
nlev_SD = [3 3 4 4 5];

% the filter used in the multiscale decomposition
smooth_func = @rcos;

% Pyr_mode can be set to be 1, 1.5 or 2. Correspondingly, the redundancy 
% factors of the transform will be 2.33, 1.59, or 1.33, respectively.
Pyr_mode = 1;

% the filter used in the directional decomposition
dfilt = 'pkva';

x = zeros(512);

%% ContourletSD Decomposition
y = ContourletSDDec(x, nlev_SD, Pyr_mode, smooth_func, dfilt);

%% Display one directional bandpass basis image
scale = 4;
d_idx = 15;  %% index for the directional subband

subband = y{scale}{d_idx};
sz = size(subband);
subband(sz(1) /2 , sz(2) / 2) = 1;
y{scale}{d_idx} = subband; %% put a unit impulse in this subband

%% ContourletSD Reconstruction
basis = ContourletSDRec(y, Pyr_mode, smooth_func, dfilt);

%% Show the basis in the frequency domain
figure
F = fftshift(fft2(basis));
F = abs(F(1:4:end, 1:4:end));
mesh(F);
view(-45,85);
title('One ContourletSD Basis Image in the Frequency Domain');

%% Show the basis in the spatial domain
figure
mesh(-basis(120:350,120:350))
view(-45,90);
title('One ContourletSD Basis Image in the Spatial Domain');