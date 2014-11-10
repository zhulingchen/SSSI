function c = joinSectionName( obj, sn )
%
%function c = joinSectionName( obj, sn )
%   returns char, eg1. '~Ascii', eg2. '~Log_data[1] | Log_definition[1]'
%   sn = cellstr array such as obj.sections{1}{1}
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

if ischar(sn)
    sn = {sn};
end

if isempty(sn) || numel(sn) >2
    warning('crewes:las:joinsectionname',...
        'Input may not be a section name, returning empty string');
    c = '';
else
   if isequal(numel(sn), 1)
       c = char(sn);
   elseif strncmpi(sn{1},'~ascii',6) || strncmpi(sn{1},'~a ',3);
       c = sn{1};       
   else
       c = sprintf('%s | %s',sn{:});
   end    
end

end %end joinSectionName()