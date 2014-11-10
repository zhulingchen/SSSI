function sigout=removenull(sigin,z,nul,inter,cons)
% sigout=removenull(sigin,z,nul,inter,cons)
%
% REMOVENULL will remove any null values from a log and replace it either
% with a nan or interpolate the missing data.  If the log has large areas of
% null values the interpolation will be inaccurate.
%
% sigin = signal or log input
%  z    = the depth or time vector coresponding to sigin
% nul   = the null value **DEFAULT =-999.25**
% inter = the type of interpolation
%            0= no interpolation the null value will be replaced with nan
%            1= linear interploation will be applied
%            2= pchip interplotation will be applied **DEFAULT**
%            3= spline interpolation will be applied
%            4= a constant will replace the null values
% cons  = a constant to fill the null values with when inter=4
%
% sigout = signal with nulls removed and interpolated depending on settings
%
% H.J.E. Lloyd November 2013
%
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.
 
% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
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
if nargin <4
    inter=2;
end
if nargin <3
nul=-999.25;
end
if inter==4
    if nargin<5
        error('Constant must be specified if interpolation method is set to 4');
    end
end
method=['linear';'pchip ';'spline'];
indnull=find(sigin==nul);
indgood=find(sigin~=nul);

sigout=sigin;
if ~isempty(indnull)
if inter==0
    sigout(indnull)=nan;
elseif inter==4
    sigout(indnull)=cons;
else
    sigout(indnull)=interp1(z(indgood),sigin(indgood),z(indnull),method(inter,:));
end
end

