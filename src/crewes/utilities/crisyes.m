function oneforyes = crisyes(txt)
% CRISYES  Tests text to see if it looks like 'Yes'.  If so, returns 1, else 0.

    if ischar(txt)
        switch lower(txt)
        case {'y','yes','true','1','affirmitive','oui','si','yup','toowoomba'}
            oneforyes = 1;
        otherwise
            oneforyes = 0;
        end
    elseif isnumeric(txt)
        oneforyes = (txt ~= 0);
    else
        error('Input to crisyes should be character text');
    end
    return

