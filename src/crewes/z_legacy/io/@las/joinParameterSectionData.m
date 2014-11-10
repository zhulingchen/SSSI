function ca = joinParameterSectionData( obj, ca )
%
% function ca = joinParameterSectionData( obj, ca )
%   returns formatted char from a cell array
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

[m, n] = size(ca);  %number of data columns
if isequal(m,1)
    % should have 6 rows for a parameter section, unless it's '~other'
    return
end

%size(ca)
% Add '.' to start of units
ca(2,:)=strcat('.', ca(2,:));

% Add '"' to start and end of anything containing obj.delimiter
% d = strfind(c(3,:),obj.delimiter)
% ** TODO **

% Add ':' to start of description
ca(4,:)=strcat({' : '},ca(4,:));

% Add {} to start of format if there is one
t = cellfun(@(X) ~isempty(X),ca(5,:));
ca(5,t)=strcat(' {', ca(5,t), '}');

% Add | to start of definition if there is one
t = cellfun(@(X) ~isempty(X),ca(6,:));
ca(6,t)=strcat({' | '}, ca(6,t));

% Get maximum number of characters in each column
nchar = cellfun(@(X) length(X), ca);
nchar = max(nchar,[],2);

% Check for minimum value size, if zero, force to 5 characters
if nchar(3) < 10
    nchar(3)=10;
end

% Check for minimum description size, if zero, force to 5
% characters
if nchar(4) < 10
    nchar(4)=10;
end

for i = 1:n
    col=ca(:,i); % for each column in the cell array
    %b{1} = mnemonic
    %b{2} = units
    %b{3} = value
    %b{4} = description
    %b{5} = format
    %b{6} = definition
    cao{i} = ...
        sprintf('%1$*7$s%2$*8$s   %3$*9$s%4$*10$s%5$*11$s%6$*12$s',...
        col{1}, col{2}, col{3}, col{4}, col{5}, col{6},...
        -nchar(1), -nchar(2), -nchar(3),...
        -nchar(4), -nchar(5), -nchar(6) );
end

ca = cao;
    
end %end function

