function dclose = closeline(xxx,yyy,xw,yw) 
% function dclose = closeline(xxx,yyy,xw,yw) 
%
% function dclose = closeline(xxx,yyy,xw,yw) 
% function linewelltie(action) 
%  finds closest possible point in line area
%    (line defined by xxx(i),yyy(i))
%     to given well location (xw,yw). 
%    Line area is defined by
%     circumscribed rectangle which surrounds line.
%  returns dclose, the distance of closest 
%     possible point in the line area to the well.
%
%  T. N. BISHOP,  DECEMBER 1993,  CPTC CANADA
%   see also seis2well, linewelltie
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
    xxxmin = min(xxx);
    yyymin = min(yyy);
    xxxmax = max(xxx);
    yyymax = max(yyy);
    if xw >= xxxmax
      xclose=xxxmax;
    elseif xw <= xxxmin
      xclose=xxxmin;
    else
      xclose=xw;
    end
    if yw >= yyymax
      yclose=yyymax;
    elseif yw <= yyymin
      yclose=yyymin;
    else
      yclose=yw;
    end
    dclose=sqrt( (xclose-xw).^2 + (yclose-yw).^2 );
