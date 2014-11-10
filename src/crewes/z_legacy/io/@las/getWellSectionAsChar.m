function c = getWellSectionAsChar(obj,sn,outputversion)
%
% function c = getWellSectionAsChar(obj,sn,outputversion)
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

if isequal(inputversion,outputversion)
    c = obj.getSectionAsChar(sn);
else
    ca = obj.getSection(sn);
    ca = swapColumns(ca,inputversion,outputversion);
    ca = updateWellSection(ca,outputversion);
    c  = obj.getSectionAsChar(ca,outputversion);
end
    
end % end function

function ca = swapColumns (ca,inputversion,outputversion)
    if inputversion < 2.0 && outputversion < 2.0
        %nothing to do
    elseif inputversion < 2.0 || outputversion < 2.0
        idx = isnan(str2double(ca{3}(3,:)));
        tmp = ca{3}(4,idx);
        ca{3}(4,idx)=ca{3}(3,idx);
        ca{3}(3,idx)=tmp;        
    end
end

function ca = updateWellSection (ca,outputversion)

%ca{1} = {'~Well'};

%set comments (if any) to null
%ca{2} = cell(1,0);  

%find minimum required mnemonics
idx = strcmpi(ca{3}(1,:),'STRT') | ...
      strcmpi(ca{3}(1,:),'STOP') | ...
      strcmpi(ca{3}(1,:),'STEP') | ...
      strcmpi(ca{3}(1,:),'NULL') | ...
      strcmpi(ca{3}(1,:),'COMP') | ...
      strcmpi(ca{3}(1,:),'WELL') | ...
      strcmpi(ca{3}(1,:),'FLD')  | ...
      strcmpi(ca{3}(1,:),'LOC')  | ...
      strcmpi(ca{3}(1,:),'PROV') | ...
      strcmpi(ca{3}(1,:),'SRVC') | ...
      strcmpi(ca{3}(1,:),'DATE') | ...
      strcmpi(ca{3}(1,:),'UWI');

%strip everything but required mnemonics
ca{3}(:,~idx)=[];

if outputversion >= 3.0 % Add required mnemonics not present in LAS < 3.0
    ca{3} = horzcat(ca{3}, ...
    {...
        'CTRY','PROV','LIC','LATI','LONG','GDAT'; ...
        ''    ,''    ,''   ,''    ,''    ,''    ; ...
        ''    ,''    ,''   ,''    ,''    ,''    ; ...
        'COUNTRY','PROVINCE','LICENSE NUMBER',    ...
            'X LOCATION','Y LOCATION','Geodetic Datum'; ...
        ''    ,''    ,''   ,''    ,''    ,''    ; ...
        ''    ,''    ,''   ,''    ,''    ,''    ;} ...        
    )
        
        
end

end %end function updateVersion
        
