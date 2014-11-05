function VelocityModel = layerModelGenerator(Vel, LayerSize, Width, varargin)
% This function genearte a layered cake velocity model
% Input:
% Vel:  velocity of each layer in unit of m/s
% Width: width of the velocity model in unit of samples
% LayerSize: depth of each layer in unit of samples
% Output:
% VelocityMoel: Velocity model in size of sum(LayerSize) by Width in which
% the velocity of each layer is specified by Vel
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lijun Zhu (gatechzhu@gmail.com) Center for Energy & Geo 
% Processing Georgia Institute of Technology

%% Input check
if length(Vel) ~= length(LayerSize)
    error('The length of Vel and that of LayerSize need to be the same');
end

if length(varargin) < 1
    Length = 1;
else
    Length = varargin{1};
end

%% Generate velocity model
% initialize
VelocityModel = zeros(sum(LayerSize), Width, Length);

% specify velocity
for i = 1 : length(Vel)
    from = sum(LayerSize(1:(i-1)))+1;
    to = sum(LayerSize(1:i));
    
    VelocityModel(from:to,:,:) = ones(LayerSize(i),Width, Length) * Vel(i);
end
return
end % EOF