function segy = SEGY_ReplaceTrace(segy, tnum, newtr)

% segy = SEGY_ReplaceTrace(segy, tracenum, newtr)
% 
% replaces the data in the trace 'tracenum' with newtr. It doesn't change
% any headers or anything like that, it just drops in the new data as a
% replacement for the old.
%
% Chad Hogan, 2008
%
% $Id: SEGY_ReplaceTrace.m,v 1.1 2008/03/04 22:38:14 cmhogan Exp $

trace = SEGY_ReadTrace(segy, tnum);

if (length(trace.data) ~= length(newtr))
    error('Your new trace must be the same size as the existing trace!');
end

trace.data = newtr;

SEGY_TraceSeek(segy, tnum);
fseek(segy.FILE, 1920/8, 0);

if(fwrite(segy.FILE, newtr, 'float') ~= length(newtr)) 
    disp(ferror(FILE)); 
    error('trace replace failed'); 
end