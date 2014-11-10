function s = getSectionAsChar(obj,sn,delimiter)
%
%function s = getSectionAsChar(obj,sn,delimiter)
%  returns a formatted character array 
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
%

if isequal(nargin,2)
    delimiter = obj.delimiter;
end

if iscell(sn) && isequal(numel(sn),4)
    ca = sn;
else
    sidx = obj.getSectionIndices(sn);
    
    if isempty(sidx)
        s = cell(1,0);
        return
    end
    
    ca = obj.getSection(sidx);
end
% section name
s1 = obj.joinSectionName(ca{1});

% section comments
s2 = horzcat(s1,ca{2});

% data section
if obj.isDataSection(ca{1})
    s3 = obj.joinDataSectionData(ca{3},delimiter);
    s = char(vertcat(s2',s3));
% parameter section    
else
    s3 = obj.joinParameterSectionData(ca{3});
    s = char(horzcat(s2{1},s3));
end

end %end function