function trimmed=nulltruncate(arraywithnulls)
% 
% NULLTRUNCATE
% 
% trimmed = nulltruncate(arraywithnulls);
%
% This function takes an array, and chops it off at the first occurence of
% a NULL (integer 0). It does NOT call DEBLANK or STRUNPAD.
%
% Chad Hogan, 2004.

% $Id: nulltruncate.m,v 1.4 2006/04/26 17:31:02 henry Exp $

nullindex = find(~arraywithnulls);
if ~isempty(nullindex)
    index = nullindex(1);
    trimmed = arraywithnulls(1:index - 1);
else
    trimmed = arraywithnulls
end

