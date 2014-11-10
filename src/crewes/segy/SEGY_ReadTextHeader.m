function thead = SEGY_ReadTextHeader(segy, format)

% thead = SEGY_ReadTextHeader(segy, format)
% 'segy' is a SEGY structure, open for reading by SEG_OpenFile().
% 'format' is either 'ascii' or 'ebcdic', and specifies the format in which
% the SEGY header was written.
% thead will be returned in ascii format.
%
% Chad Hogan, 2008
%
% $Id: SEGY_ReadTextHeader.m,v 1.2 2009/07/24 15:49:19 kwhall Exp $


if nargin < 2
    format = 'ascii';
end

% Make sure we're at the start of the file
% Assume only one textual header for now
if(fseek(segy.FILE, 0, 'bof') ~= 0)
    % rewind the file to the textual header start
    disp(ferror(FILE));
    error('seeking failed');
end

% First we'll grab a template header to fill in
thead = SEGY_GetTextHeader;

thead = strcat(fread(segy.FILE, 3200, 'char'));

if strcmp(lower(format), 'ebcdic')
    thead = ebcdic2ascii(thead);
end