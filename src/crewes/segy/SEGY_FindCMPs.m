function segy = SEGY_FindCMPs(segy, binsize)

% segy = SEGY_FindCMPs(segy, binsize)
%
% Zips through a segy file and figures out how many unique CMPs there are
% given a 'binsize' to choose how loose the exact CMP location can be.
% This code assumes a LINEAR SURVEY. If your survey is deviated, you're
% going to have to write your own CMP finder, probably using SEGY_FindIndex.
% Thems the breaks.
%
% cmp distances (linear, starting at the beginning of the survey) are
% returned in segy.cmps
%
% Chad Hogan, 2008
%
% $Id: SEGY_FindCMPs.m,v 1.1 2008/03/04 22:36:27 cmhogan Exp $
bin = binsize/2;

disp(['Going to search ' num2str(segy.numtraces) ' traces']);

cmps = [];

for idx = 1:segy.numtraces
    SEGY_TraceSeek(segy, idx);
    fseek(segy.FILE, 72, 0);    % move to sx
    thissx = fread(segy.FILE, 1, 'int');
    thissy = fread(segy.FILE, 1, 'int');
    thisgx = fread(segy.FILE, 1, 'int');
    thisgy = fread(segy.FILE, 1, 'int');

    cmpx = (thissx + thisgx) / 2;
    cmpy = (thissy + thisgy) / 2;
    
    if (idx == 1)
        x0 = cmpx;
        y0 = cmpy;
    end
    
    dx = cmpx - x0;
    dy = cmpy - y0;
    
    dist = sqrt(dx^2 + dy^2);
    
    dist = round(dist / bin) * bin;
    
    if(~length(find(cmps == dist)))
        cmps(end+1) = dist;
    end
    
    if(mod(idx, 1000) == 0)
        disp(['done ' num2str(idx) ' of ' num2str(segy.numtraces)]);
    end
end

segy.cmps = cmps;

disp(['Found ' num2str(length(segy.cmps)) ' CMP locations']);