function format = guessTextFormat(text)
%function format = guessTextFormat(text)
% Guess if text header is ascii or ebcdic format
% Theory:
%   EBCDIC characters can have values > 127, Standard ASCII cannot
% Returns:
%   'ascii' unless it can be proven to be 'ebcdic'

format = 'ascii'; %assume ascii

if any(text > 127)
    format = 'ebcdic';
end

end % end guessTextFormat

