function [h, t, gather] = SEGY_ReadCMPGather(segy, cmpnum, binsize)

% [h, t, gather] = SEGY_ReadCMPGather(segy, cmpnum, binsize)
%
% Returns the CMP gather binned with binsize at segy.cmps(cmpnum)
%
% Chad Hogan, 2008
%
% $Id: SEGY_ReadCMPGather.m,v 1.1 2008/03/04 22:37:12 cmhogan Exp $

if isnan(segy.cmps)
    error('You did not find CMPs! Use SEGY_FindCMPs() to find them.');
end

if (cmpnum > length(segy.cmps))
    error('Aint that many CMPs gathers, hotshot.');
end
bin = binsize/2;

disp(['Going to search ' num2str(segy.numtraces) ' traces']);

gsize = 0;

cmp = segy.cmps(cmpnum);

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
    
    if(dist == cmp)
        trace = SEGY_ReadTrace(segy, idx);
        gather(:, gsize + 1) = trace.data;
        thish = sqrt((thisgx - thissx).^2 + (thisgy - thissy).^2);
        if(thisgx == thissx)
            thish = sign(thisgy - thissy) * thish;
        else
            thish = sign(thisgx - thissx) * thish;
        end
        h(gsize+1) = thish;
        gsize = gsize + 1;
    end
    
    if(mod(idx, 1000) == 0)
        disp(['done ' num2str(idx) ' of ' num2str(segy.numtraces)]);
    end
end

t = ((1:size(gather, 1))-1) *segy.bhead.hdt;
