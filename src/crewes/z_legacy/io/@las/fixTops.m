function obj = fixTops(obj)
%
% function obj = fixTops(obj)
%
% Section ~Tops is not formally defined in the LAS standard for 
% versions prior to 3.0, so different companies did different things
%
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

if inputversion < 3.0 && sum(obj.isDataSection('~t'))
    % Tops Section was not recognized as a parameter section on input,
    % need to fix. Possible to avoid this in splitLas() ?
    sn = obj.getSectionNames('~t');
    
    for m = 1:numel(sn) %take care of multiple ~t sections LAS ver < 3
        if obj.isDataSection(sn{m})
            ca = obj.getSectionData(sn{m});
            
            % which cells _do_ convert to doubles ?
            k = cellfun(@(x) ~isnan(str2double(x)),ca);
            
            % set numbers to empty strings, and join
            tn = ca;
            tn(k) = cellstr(blanks(0));
            tn = strtrim(obj.joinDataSectionData(tn,' '));
            
            %sheesh, never ends. remove multiple spaces from within log mnemonics
            % re-write in joinDataSection?
            tn = regexprep(tn,'\s+',' ');
            
            % get top depths (td)
            td = ca(:,1); %set td to the first column of ca
            
            for n = 1:size(ca,2) % column by column
                k = cellfun(@(x) ~isnan(str2double(x)),ca(:,n));
                td(k) = ca(k,n); %accumulate numerics in td
            end
            
            %create new data cell array filled with empty strings
            tdat = cell(6,length(tn));
            tdat(:) = cellstr(blanks(0));
            
            % set top names (tn) and top depths (td)
            tdat(1,:) = tn;
            tdat(3,:) = td;
            
            % reset tops section data in object
            obj.sections = {'sectiondata',sn{m},tdat};  %update section data
            obj.sections = {'sectiontype',sn{m},false}; %set parameter section
        end
    end %end for
end

end %end function