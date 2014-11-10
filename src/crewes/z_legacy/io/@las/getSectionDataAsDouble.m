function v = getSectionDataAsDouble(obj,sn)
%
%function ca = getSectionDataAsDouble( obj, sn )
% Returns cell array(s) of numeric data selected from obj.sections based on 
% indices contained in sn (if numeric) or on case-invariant matches between 
% section names in obj.sections and the contents of sn (char or cellstr).
%
% Note that this function can return the contents of more than one section:
% eg. sn='~log' could match '~log_parameter','~log_definition','~log_data',
% etc.
%
% If no matches to sn are found in the object, getSectionData returns an empty
% cell array
%
% If the cell array contains char() data, getSectionDataAsDouble replaces
% it with str2double(obj.lognull)
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

v = obj.getSectionData(sn);
if isempty(v)
    v=cell(1,0);
    return
end

if iscellstr(v)
    v = str2double(v);
    v(isnan(v)) = str2double(obj.lognull);
    v = {v};
elseif iscell(v)
    for idx = 1:numel(v)
        t = str2double(v{idx});
        t(isnan(t)) = str2double(obj.lognull);
        v{idx} = t;        
    end
end

end


