function texthdr = SEGY_ModifyTextHeaderLine(line,comment,texthdr)
    % Function to modify row of text header provided by
    % texthdr = SEGY_GetTextHeader()
    %
    % line  = text header row to modify (1-40)
    % comment = string to embed in text header
    % thline  = modified string to embed in text header
    % texthdr = header to modify
    
    if(line > 40)
        error('Text Header cannot have more than 40 rows');
    end
    i = (line-1)*80+1; %line*ncols +1
    thline = sprintf ('C%2d %s',line,comment);
    if (length(thline) > 80)
        warning('Modified comment is too long! Truncating...');
        thline = thline(1:80);
    end    
    
    texthdr(i:i+length(thline)-1) = thline;
end
