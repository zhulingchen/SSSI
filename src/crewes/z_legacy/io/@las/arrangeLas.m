function obj = arrangeLas(obj)
%
% function obj = arrangeLas(obj)
% Arrange regexp results from obj.splitlas into a cell array of sections
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

% grab all section names from obj
s = obj.sections(1,:);

% find indices of non-empty section names in cell array
idx = find(cellfun(@(x) ~isempty(x), s)==1);

% separate any associtated section names based on '|' (LAS 3)
sn = regexp(s(idx),'\s\|\s','split');

% find indices for start and end of each section in cell array
startidx = idx+1;  %section starts one row after section name
endidx = idx-1;    %section ends one row before section name
endidx(end+1) = length(s); %except for last section, which ends at EOF
endidx = endidx(2:end); %and starts at index 1, not 0

%pre-allocate arrays
fc{length(idx)}=[];

np(idx)=0;
nd(idx)=0;

for i = 1:length(idx) % i = sequential section number within las file
    %%set section name
    fc{i}{1} = sn{i};
    
    %get section contents
    sc = obj.sections(:,startidx(i):endidx(i)); %section contents

    % How many parameters or data rows? mnemonics should be in row 3, data
    % in row 4
    np(i) = sum(cellfun(@(X) ~isempty(X),sc(3,:)));
    nd(i) = sum(cellfun(@(X) ~isempty(X),sc(4,:)));
    
    if isempty(sc)
        %section has no information in it
        fc{i}{2}={[]};
        fc{i}{3}={[]};
    else
        % get comment lines
        comments = sc(2,:);
        c=cellfun(@(X) isempty(X),comments);
        comments(c) = [];
        fc{i}{2}=comments;
        %nuke empty cols corresponding to comment lines in sc
        sc(:,~c) = [];
        
        if ~isequal(np(i),0) %no data lines in this section
            disp(['     Section: ',sn{i}{1}])
                %'parameter section', remove section name and data rows;
                sc(1:2,:)=[]; %remove section name and comment lines
                sc(2,:)=[]; %remove comment lines
                
                %find exponential formats
                idx = cellfun(@(X) ~isempty(X),regexp(sc(5,:),'E'));
                
                %fix exponential formats (eg. E00.00E+00 => E0.000E+00)
                sc(5,idx) = cellfun(@(X) obj.modifyExponentialFormat(X),...
                    sc(5,idx),'UniformOutput',false);
                
                fc{i}{3} = sc;
                fc{i}{4} = false; % not data section

        else
            
                %Three cases, '~other', '~ascii',
                % '~yyy_data' | yyy_definition'
                if strncmpi(sn{i}{1},'~o',2)
                    disp(['     Section: ',sn{i}{1}])
                    fc{i}{3} = sc(4,:);
                    fc{i}{4} = false; %not data section
                else
                    disp(['     Section: ',sn{i}{1}])                 
                    fc{i}{3} = sc(4,:);
                    fc{i}{4} = true; %data section
                end                
        end
    end
end
obj.sections = {'all',fc};
obj.sectionNames = obj.setAssociatedSections(sn);

end

