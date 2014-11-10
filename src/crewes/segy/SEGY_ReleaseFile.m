
function status = SEGY_ReleaseFile(segy)

% SEGY_ReleaseFile(segy)
%
% cleans up the opened file and so on. Just run this before you clear or
% something. I don't think it's terribly crucial but I suppose it's
% considerate.
%
% Chad Hogan, 2008
%
% $Id: SEGY_ReleaseFile.m,v 1.2 2009/07/24 15:49:19 kwhall Exp $


fclose(segy.FILE);
