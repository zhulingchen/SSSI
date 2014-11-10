function v = getMnemonicFormat(obj,sn,mn)
%
%function v = getMnemonicFormat(obj,sn,mn)
% sn = (partial) section name (char)
% mn = LAS mnemomic (char)
% v  = Mnemonic description (char)
%
% mnemonics are in the 1st row
% units are in the 2nd row
% values are in the 3rd row
% descriptions are in the 4th row
% formats are in the 5th row
% associations are in the 6th row
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

% Get indices for mnemonic of interest in the section data
midx = obj.getMnemonicIndex(sn,mn);

switch sum(midx)
    case 0 %mnemonic not found in first row of cellarray
        v='';
    otherwise %formats are in the 5th row
        ca = obj.getSectionData(sn);        
        v=ca{5,midx};
end

end
