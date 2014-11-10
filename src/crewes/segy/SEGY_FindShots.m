function segy = SEGY_FindShots(segy)

% segy = SEGY_FindShots(segy)
%
% Finds unique shot locations segy.sx and segy.sy. Also sets the number of
% traces for that shot in segy.shottraces.
%
% Chad Hogan, 2008
%
% $Id: SEGY_FindShots.m,v 1.1 2008/03/04 22:36:42 cmhogan Exp $

disp(['Going to search ' num2str(segy.numtraces) ' traces']);

sx = [];
sy = [];
shottraces = [];

nt = 0;

for idx = 1:segy.numtraces
    SEGY_TraceSeek(segy, idx);
    fseek(segy.FILE, 72, 0);    % move to sx
    thissx = fread(segy.FILE, 1, 'int');
    thissy = fread(segy.FILE, 1, 'int');
    nt = nt + 1;

    if (idx == 1)
        lastx = thissx;
        lasty = thissy;
    end
    
    % If this is a unique shot that we haven't found before
    if(lastx ~= thissx || lasty ~= thissy)
        lastx = thissx;
        lasty = thissy;
        sx(end+1) = thissx;
        sy(end+1) = thissy;
        shottraces(end+1) = nt;
        nt = 0;
    end
    if(mod(idx, 1000) == 0)
        disp(['done ' num2str(idx) ' of ' num2str(segy.numtraces)]);
    end
end
disp(['Found ' num2str(length(sx)) ' unique shot locations']);

segy.sx = sx;
segy.sy = sy;
segy.shottraces = shottraces;