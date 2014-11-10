function obj = read(obj)
%function obj = read(obj)
%  obj.fileName and obj.fullFileName should exist
%  functions reads all text lines from LAS file into obj.sections
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

obj.sections=[];

try
    disp(['Reading file: ' obj.fileName])
    fid = fopen(obj.fullFileName,'r');
    
    % Read every line of LAS file into memory
    s=textscan(fid,'%s','delimiter','\n');
    
    fclose(fid);
    
    % simplify cell array structure
    s = [s{:}];
    
    % remove empty rows
    s(cellfun(@isempty,s))=[];
    
    obj.sections={'all',s};
catch ex
    warning('crewes:las:read',ex.message);
end
    
  
end