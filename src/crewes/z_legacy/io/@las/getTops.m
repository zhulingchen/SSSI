function [tn, td] = getTops(obj,tparm,tdef,tdat)
%
% function [tn, td] = getTops(obj,tp,tdef,tdat)
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

%assume failure
tn='';
td=[];

%inputversion < 3
if inputversion < 3.0
    tops_parameter = tparm;
    tops_dat = obj.getSectionData(tops_parameter);
    if isempty(tops_dat)
        return
    end
    tn = char(tops_dat(1,:));
    td = str2double(tops_dat(3,:))';
elseif inputversion >= 3.0
    %         tops_parameter  = '~tops_parameter';
    %         tops_definition = '~tops_definition';
    %         tops_data       = '~tops_data';
    
    tops_dat = obj.getSection(tdat);
    if isempty(tops_dat)
        return
    end
    %get tops names
    tn  = char(tops_dat{3}(:,...
        obj.getMnemonicIndex(tdef,'topn')));
    %get tops
    td  = str2double(tops_dat{3}(:,...
        obj.getMnemonicIndex(tdef,'topt')));
    
    %convert tops elevations to depths
    if strncmpi(obj.getMnemonicValue(tparm,'topdr'),...
            'Subsea',6)
        %get reference elevation
        % THIS needs work; check DREF for starters...
        eref = obj.getMnemonicValue('~parameter','eref');
        
        if isempty(eref)
            eref = obj.getMnemonicValue('~log_parameter','eref');
        end
        
        %convert subsea to depth
        eref = str2double(eref);
        
        if ~isnan(eref)
            td = eref -td;
        end
    end
    
    
    
end
        
end %end function
