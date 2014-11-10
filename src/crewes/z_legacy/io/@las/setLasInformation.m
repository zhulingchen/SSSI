function obj = setLasInformation( obj )
%
%function obj = setLasInformation( obj )
%
% Set basic information in obj
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

%Get information from ~version section
sidx = obj.getSectionIndices('~v');
obj.version   = obj.getMnemonicValue(sidx,'vers');
obj.delimiter = obj.getMnemonicValue(sidx,'dlm');
%obj.wrapped   = obj.getMnemonicValue(sidx,'wrap');
obj.wrapped = 'NO';
obj.sections  = {'mnemonicvalue',sidx,'wrap','NO'};

%Get information from ~well section
obj.lognull = obj.getMnemonicValue('~w','null');

end

