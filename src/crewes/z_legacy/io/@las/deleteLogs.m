function obj = deleteLogs( obj, lm, sn )
% 
% function obj = uideletelogs(obj,lm,sn)
%
% DELETELOGS is a utility designed to allow the deletion of logs
% in LAS files. This is useful for files that contain large numbers of logs
% because LOGEDIT does not work well with more than about 20 logs.
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

%Check that lm is a vector
if ~isvector(lm)
    warning('crewes:las:deletelogs',...
            'Input lm must be a vector; giving up');
    return    
end

%Get section indices from object
sidx = obj.getSectionIndices(sn);

%Make sure sn only matches one section name
if(numel(sidx)>1)
    warning('crewes:las:deletelogs',...
        ['Refusing to delete logs from multiple sections; giving up']);
    return
end

%Make sure the section name matched is for a data section
if ~obj.isDataSection(sidx)
    if isnumeric(sn)
    warning('crewes:las:deletelogs',...
        ['Section ' int2str(sn) ' is not a data section; giving up']);        
    else
    warning('crewes:las:deletelogs',...
        ['Section ''' sn ''' is not a data section; giving up']);
    end
    
    return
end

%Get section indices for the data section and it's associated definition
%section in the las obj, and a list of mnemonics from the definition section

asidx  = cell2mat(obj.sectionNames(4,sidx));
ascamn = obj.getParameterSectionMnemonics(asidx);

%make sure we're dealing with a cellstr of lognames
if islogical(lm) %input lm is logical
    if numel(lm) ~= numel(ascamn)
        warning('crewes:las:deletelogs',...
            ['Logical indices in lm specify a different number of logs '...
             'than contained in sn, giving up']);
    else
        b = lm;
    end
elseif ischar(lm) || iscellstr(lm)
    if ischar(lm)
        lm = cellstr(lm);
    end
    
    %See if the input lm match any of the mnemonics (case invariant)
    b = false(1,numel(ascamn)); %By default,don't delete any logs
    for n =1:numel(lm)
        b = b | strcmpi(lm{n},ascamn); %Set to true if lm{n} is found in asca
    end
elseif isnumeric(lm)
    %check for limits
    if min(lm) < 1 || max(lm) > numel(ascamn)
        warning('crewes:las:deletelogs',...
            'Indices provided in lm are out of range, giving up');
        return
    end
    b = false(1,numel(ascamn)); %By default,don't delete any logs
    b(lm) = true;    
else
    warning('crewes:las:deletelogs',...
        ['Input log names must be of type logical, numeric, char or '...
         ' cellstr, giving up']);
    return
end

%At this point, any input should have been reduced to a logical vector

if sum(b) > 0 %We have some logs to delete!
    %get section data from object
    sn_dat = obj.getSectionData(sidx);
    asn_dat = obj.getSectionData(asidx);
    
    %Make sure we aren't deleting depth or time axes
%     b = strncmpi(logs_table(:,2),'dept',4);
%     b = b | strncmpi(logs_table(:,2),'etim',4);
    
    %set indices matching true in 'b' to []
    sn_dat(:,b)=[];
    asn_dat(:,b)=[];
    
    %update section data in object
    obj.sections={'sectiondata',sidx,sn_dat};
    obj.sections={'sectiondata',asidx,asn_dat};
else
    warning('crewes:las:deletelogs','Did not delete any logs');
end

end %end function



