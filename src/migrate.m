function m = migrate(travelTime, xRecGrid, xShotAndRecGrid, shot, dt, nz, ixs)
% Migrate a shot record for a given travel time between shot (source) and
% gather (receiver) using a simple Kirchoff Migration algorithm
%
% Inputs:
%   travelTime      travel time array
%   shot            shot array (nz,nx)
%   dt              sampling time
%   nz              number of samples in z direction
%   ixs             shot location in shot
%
% Outputs:
%   m               migrated image (nz,nx)

% Copyright 2010 The MathWorks, Inc.
% All rights reserved

nx = size(shot,2);
m = zeros(nz,nx);

for ixr = 1:length(xRecGrid)
    xr = xRecGrid(ixr);
    [~, idxXRec] = ismember(xr, xShotAndRecGrid);
    
    it = shot2RecTime(travelTime, ixs, idxXRec, dt, nx);
    m = m + reshape( shot(it, xr), nz, nx);
end