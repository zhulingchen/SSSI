function t = isLas(obj)
%
%function t = islas(obj)
% Scans 'filename' for the first line that is not a comment,
% then checks if that line starts with '~Version'. If it does
% islas() returns true.
%
%  'filename' = character string
%  't'        = true or false

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geoscience of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by 
% its author (identified above) and the CREWES Project.  The CREWES 
% project may be contacted via email at:  crewesinfo@crewes.org
% 
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) Use of this SOFTWARE by any for-profit commercial organization is
%    expressly forbidden unless said organization is a CREWES Project
%    Sponsor.
%
% 2) A CREWES Project sponsor may use this SOFTWARE under the terms of the 
%    CREWES Project Sponsorship agreement.
%
% 3) A student or employee of a non-profit educational institution may 
%    use this SOFTWARE subject to the following terms and conditions:
%    - this SOFTWARE is for teaching or research purposes only.
%    - this SOFTWARE may be distributed to other students or researchers 
%      provided that these license terms are included.
%    - reselling the SOFTWARE, or including it or any portion of it, in any
%      software that will be resold is expressly forbidden.
%    - transfering the SOFTWARE in any form to a commercial firm or any 
%      other for-profit organization is expressly forbidden.
%
% END TERMS OF USE LICENSE

t=false;

fstats = dir(obj.fullFileName);

if isempty(fstats)
    return
end

if isequal(fstats.bytes,0)
    warning('crewes:las:islas',...
        'file is empty');
    return
end
   
try
    fid = fopen(obj.fullFileName,'r'); %open filename for reading
    
    % Read the filename line by line into memory, skipping comment
    % lines until we find a non-empty line
    row = textscan(fid,'%s',1,'delimiter','\n','CommentStyle','#');
    while isempty(char(row{1}))
        row = textscan(fid,'%s',1,'delimiter','\n','CommentStyle','#');
    end
    
    fclose(fid);
    
    % Check if line starts with '~Version'
    %If it does not, this is not a valid LAS filename
    v = '~v';
    t = strncmpi(strtrim([row{:}]),v,length(v));
    
    if ~t
        warning('crewes:las:islas', ...
            'file is not LAS format');
    end
catch ex
    warning('crewes:las:islas',ex.message);
    t = false;
end %end try

end %end function