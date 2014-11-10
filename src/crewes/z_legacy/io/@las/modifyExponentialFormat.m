function fmt = modifyExponentialFormat (fmt)
%
%function fmt = modifyExponentialFormat (fmt)
%  modifies LAS exponential format so that, for example:
%     'E000.0E+00' becomes 'E0.000E+00'
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

if isempty(regexp(fmt,'E','ONCE'))
    warning('crewes:las:modifyexponentialformat',...
        [fmt ' does not appear to be a LAS exponential format string']);
    return
end

% Find index number for decimal point in the format string
idx = regexp(fmt,'\.','start');

% if the format is 'correct', ie E0.0.... the decimal point is at index 3
if ~isequal(idx,3)
    % but that would be too much to hope for; change it.
    fmt(idx)='0';
    fmt(3)='.';
end
    
end