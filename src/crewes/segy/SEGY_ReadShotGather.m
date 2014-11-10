function [x, t, shot] = SEGY_ReadShotGather(segy, snum)

% [x, t, shot] = SEGY_ReadShotGather(segy, snum)
%
% Reads a shot gather from 'segy' dataset (as opened by SEGY_OpenFile). 
% Returns the selected shotgather from shot number 'snum'. 'snum' is just
% a sequential number. If there are 200 unique shot locations, it'll return
% the data from the 'snum'th shot. 
%
% Chad Hogan, 2008
%
% $Id: SEGY_ReadShotGather.m,v 1.1 2008/03/04 22:37:36 cmhogan Exp $

if isnan(segy.sx)
    error('You did not find shots! Use SEGY_FindShots() to find them.');
end

if (snum > length(segy.sx))
    error('Aint that many shots, hotshot.');
end

skiptraces = 0;

if(snum > 1)
   skiptraces = sum(segy.shottraces(1:(snum-1))); 
end

shot = zeros(segy.bhead.hns, segy.shottraces(snum));

for idx = 1:segy.shottraces(snum)
    trace = SEGY_ReadTrace(segy, skiptraces + idx);
    shot(:,idx) = trace.data;
    x(idx) = trace.gx;
    y(idx) = trace.gy;
end

t = ((1:length(trace.data))-1)*segy.bhead.hdt;
d = sqrt(x.^2 + y.^2) - sqrt(x(1)^2 + y(1)^2);
[x, idx] = sort(d);
shot = shot(:, idx);
