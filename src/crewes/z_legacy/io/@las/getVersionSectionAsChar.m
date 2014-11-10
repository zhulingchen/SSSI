function c = getVersionSectionAsChar(obj,sn,outputversion)
%
% function c = getVersionSectionAsChar(obj,outputversion)
%
%   Return formated char from cell array appropriate for a given LAS
%   version
%

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

    %assume failure
    c=[];

    inputversion = str2double(obj.version);
%     outputversion;

    if outputversion < 2.0 %eg LAS 1.2
        version     = '1.2';
        description = 'CWLS LOG ASCII STANDARD - VERSION 1.2';
    elseif outputversion >= 3.0 %eg LAS 3.0
        version     = '3.0';
        description = 'CWLS LOG ASCII STANDARD - VERSION 3.0';
    else %LAS 2
        version     = '2.0';
        description = 'CWLS LOG ASCII STANDARD - VERSION 2.0';
    end
    
    if isequal(inputversion,outputversion)
        c = obj.getSectionAsChar(sn);
    else
        ca = obj.getSection(sn);
        ca = updateVersionSection(ca,version,description);
        %obj.sections={'section',sn,ca};
        c  = obj.getSectionAsChar(ca,outputversion);       
    end
    
end % end function
    
function ca = updateVersionSection (ca,version,description)
% reset version section name
% ca{1} = {'~Version'};

%set comments (if any) to null
% ca{2} = cell(1,0);  

%update version number in VERS mnemonic
ca{3}(1,:);
vers_idx = strcmpi(ca{3}(1,:),'VERS');
ca{3}(3,vers_idx) = {version};
ca{3}(4,vers_idx) = {description};

%strip everything but mnemonics VERS and WRAP
idx = strcmpi(ca{3}(1,:),'VERS') | strcmpi(ca{3}(1,:),'WRAP');
ca{3}(:,~idx)=[];

if str2double(version) >= 3.0 % Add DLM
    idx = strcmpi(ca{3}(1,:),'VERS') | strcmpi(ca{3}(1,:),'WRAP');
    ca{3}(:,~idx)=[];
    ca{3}(:,3) = { ...
        'DLM',     ...
        '',        ...
        'SPACE',   ...
        'DELIMITING CHARACTER (SPACE TAB OR COMMA)',...
        '',        ...
        ''};
end

end %end function updateVersion
        
