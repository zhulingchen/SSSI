function c = getLogDefinitionSectionAsChar(obj,sn,outputversion)
%
% function c = getLogParameterSectionAsChar(obj,sn,outputversion)
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

    inputversion = str2double(obj.version);
    
    %assume failure
    c=[];
      
    ca = obj.getSection(sn);
%     celldisp(ca)
    ca = updateLogDefinitionSection(ca,inputversion,outputversion);
    c  = obj.getSectionAsChar(ca);
    
%     % if object version and output version are the same:
%     if isequal(inputversion,outputversion) || ...
%             (inputversion < 3.0 && outputversion < 3.0)
%         
%         disp('No parameters conversion')
%            c=obj.getSectionAsChar(log_definition);
%         
%     elseif inputversion < 3.0 && outputversion >= 3.0
%         disp ('STUB: convert from v2 to v3')
%                 
%     elseif inputversion >= 3.0 && outputversion < 3.0
%         disp ('STUB: convert from v3 to v2')      
%       
%     end

end %end function

function ca = updateLogDefinitionSection(ca,inputversion,outputversion)

%reset section name
ca{1} = {'~CURVE'};

%set comments (if any) to null
%ca{2} = cell(1,0);  

%clear log data formats
ca{3}(5,:) = {''};

%clear log mnemonic associations 
ca{3}(6,:) = {''};

%strip everything but required mnemonics
%ca{3}(:,~idx)=[];

if outputversion >= 3.0 % Add required mnemonics not present in LAS < 3.0
%     ca{3} = horzcat(ca{3}, ...
%     {...
%         'CTRY','PROV','LIC','LATI','LONG','GDAT'; ...
%         ''    ,''    ,''   ,''    ,''    ,''    ; ...
%         ''    ,''    ,''   ,''    ,''    ,''    ; ...
%         'COUNTRY','PROVINCE','LICENSE NUMBER',    ...
%             'X LOCATION','Y LOCATION','Geodetic Datum'; ...
%         ''    ,''    ,''   ,''    ,''    ,''    ; ...
%         ''    ,''    ,''   ,''    ,''    ,''    ;} ...        
%     )               
end

end