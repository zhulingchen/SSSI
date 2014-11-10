function SEGY_TraceSeek(segy, tracenum)

% SEGY_TraceSeek(segy, tracenum)
%
% Seeks to a given trace number 'tracenum'. 'segy' is a struct as returned by
% SEGY_OpenFile().
%
% Chad Hogan, 2008
%
% $Id: SEGY_TraceSeek.m,v 1.1 2008/03/04 22:38:20 cmhogan Exp $

if (tracenum > segy.numtraces)
    error('tracenum is larger than the total number of traces in the SEGY!');
end
bytelen = (1920 + segy.bhead.hns * 32) / 8;

fseek(segy.FILE, 3200 + 400 + (tracenum-1)*bytelen, -1);