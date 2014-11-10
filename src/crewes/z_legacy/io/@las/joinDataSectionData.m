function dat = joinDataSectionData( obj, dat, delimiter )
%
%function cs = joinDataSectionData( obj, dat, delimiter )
%
% Assumption: section data is already cellstr; needs to be more general.
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

if nargin < 3
    delimiter = obj.delimiter;
end

dat = obj.addDoubleQuotes(dat,delimiter);
if strcmp(delimiter,' ')
    dat = obj.padcellstr(dat);
end
% interleave columns
dat=dat';
[m, n] = size(dat);

% create 2D matrix of delimiters the same size as cs
delimiter = repmat({delimiter},m,n);

% interleave 'ca' and 'delimiter' columns
dat = reshape([dat(:) delimiter(:)]', 2*m, [])';

% drop rightmost column (extraneous delimiters)
dat(:,2*m)=[];

% concatenate cellstr's so we have one cellstr per row.
dat = num2cell(dat,1);
dat = strcat(dat{:});

end

