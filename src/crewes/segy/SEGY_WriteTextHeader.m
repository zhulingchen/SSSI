function SEGY_WriteTextHeader(FILE, thead, format)
% SEGY_WRITETEXTHEADER(FILE, thead, format)
%
% This function writes a text header 'thead' to the file FILE with the
% 'format', where 'format' is either 'ebcdic' or 'ascii'. 'format' is an
% optional argument, and if unspecified will default to 'ascii'. You can
% retrieve a valid text header from the function GETTEXTHEADER. It is
% wise to get a text header from this function and then edit it as you
% see fit.
%
% This function will OVERWRITE whatever text header may be existing in
% the file. 
%
% Chad Hogan, 2004
% 
% $Id: SEGY_WriteTextHeader.m,v 1.1 2004/06/18 21:24:29 cmhogan Exp $

if nargin < 3
    format = 'ascii';
end

if strcmp(lower(format), 'ebcdic')
    thead = ascii2ebcdic(thead);
end
    

if length(thead) ~= 3200
    error('input text header must be exactly 3200 bytes long');
    return
end

fseek(FILE, 0, 'bof'); % rewind the file to the very start

if(fwrite(FILE, thead) ~= 3200)
    disp(ferror(FILE))
    error('writing failed');
end





