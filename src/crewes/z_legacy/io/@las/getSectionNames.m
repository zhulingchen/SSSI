function ca = getSectionNames( obj, sn )
%
%function ca = getSectionName( obj, sn )
% Returns cell array(s) selected from obj.sections based on indices
% contained in sn (if numeric) or on case-invariant matches between 
% section names in obj.sections and the contents of sn (char or cellstr).
%
% Note that this function can return the contents of more than one section:
% eg. sn='~log' could match '~log_parameter','~log_definition','~log_data',
% etc.
%
% If no matches to sn are found in the object, getSectionName returns an empty
% cell array
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

c = obj.getSectionIndices(sn);

if isempty(c)
    ca = cell(1,0);
else
    ca = obj.sectionNames(1,c);
end


end

