function segy = SEGY_ReplaceShotGather(segy, snum, data)

% segy = SEGY_ReplaceShotGather(segy, shotnum, data)
%
% Replaces the shot gather 'shotnum' with 'data'. Does not modify any
% headers or anything else -- it just replaces the trace data itself. You
% must have run SEGY_FindShots() before this will work.
%
% Note that your new shot must be the same size as the old shot. This
% *replaces* a shot in-place.
%
% Chad Hogan, 2008
%
% $Id: SEGY_ReplaceShotGather.m,v 1.1 2008/03/04 22:38:06 cmhogan Exp $


if isnan(segy.sx)
    error('You did not find shots! Use SEGY_FindShots() to find them.');
end

if (snum > length(segy.sx))
    error('Aint that many shots, hotshot.');
end

numtraces = segy.shottraces(snum);

if (size(data, 2) ~= numtraces)
    error('The new shot gather must have the same number of traces as the old one!');
end

skiptraces = 0;

if(snum > 1)
   skiptraces = sum(segy.shottraces(1:(snum-1))); 
end

for idx = 1:segy.shottraces(snum)
   segy = SEGY_ReplaceTrace(segy, skiptraces + idx, data(:, idx)); 
end