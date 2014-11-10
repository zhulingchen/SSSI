function segy = SEGY_OpenFile(fname, mode, endian, textfmt, iwarn)

% segy = SEGY_OpenFile(fname, mode, endianness, textfmt)
% 
% Returns a structure that has a filehandle, along with some useful
% information like:
%
% segy.FILE ........ the filehandle
% segy.numtraces ... the number of traces in the file
% segy.thead ....... the text header
% segy.bhead ....... the binary header
%
% This structure is ready for immediate trace extraction, using the
% SEGY_ReadTrace() function.
%
% See fopen() for details on 'mode' and 'endianness' choices.
% txtfmt = 'ascii' or 'ebcdic'
%
% The returned structure also contains some empty arrays like segy.idxs,
% and some NaN-valued elements like segy.sx, segy.sy. These are used in
% in other functions.
%
% ALL FILES OPENED WITH SEGY_OPENFILE SHOULD BE RELEASED WITH
% SEGY_RELEASEFILE
%
% Chad Hogan, 2008
%
% $Id: SEGY_OpenFile.m,v 1.3 2009/07/24 15:49:19 kwhall Exp $

if(nargin<5)
    iwarn=1;
end
if nargin < 4
    textfmt = 'ascii';
end

segy.sx = NaN;
segy.sy = NaN;
segy.cmps = NaN;
segy.shotgathersize = NaN;
segy.idxs = [];

segy.FILE = fopen(fname, mode, endian);

if(segy.FILE == -1)
    error('OPEN FAILED'); 
end

segy.thead = SEGY_ReadTextHeader(segy, textfmt);
segy.bhead = SEGY_ReadBinaryHeader(segy);

fseek(segy.FILE, 0, 'eof');
allbytes = ftell(segy.FILE);

tracebytes = allbytes - 3600;
segy.numtraces = tracebytes / segy.bhead.tracebytelen;

if(mod(segy.numtraces, 1) ~= 0 && iwarn==1)
    warning('Number of traces is not an integer! Verify that SEGY data contains proper headers etc.');
end

if(iwarn==1)
    disp(['SEGY opened. ' num2str(segy.numtraces) ' traces found, sampling at ' ...
        num2str(segy.bhead.hdt/1000) ' ms']);
end

