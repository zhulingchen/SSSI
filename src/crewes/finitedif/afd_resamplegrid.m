function [section,x,z,dxz] = afd_resamplegrid( oldsection,oldx,oldz,newdxz)
%AFD_RESAMPLEGRID is a program to convert sections that do not have equal
%gridspacing in both the x and z directions, ie. dx is not equal to dz.
%this program also will start both x and z at the zero as this is required
%for the finite differencing toolbox.  If you would like to have these
%vectors start at a different location simply add the starting point to the
%vector ie. z=z+400;
%
%[section,x,z,dxz] = afd_resamplegrid( oldsection,oldx,oldz,newdxz)
%
%   The variables that are required for this funtion are:
%        oldsection... this is the section that you would like to have
%                      resampled.
%        oldx... This is the inline vector of the section you would like to
%                have resampled.  dx will be calculated from this vector.
%        oldz... This is the depth vector of the section you would like to
%                have resampled.  dz will be calculated from this vector.
%        newdxz... This is the sampling rate that you would like in both 
%                  the depth and the inline directions.  ie. dz=dx=newdxz
%                  ***********default=10***********
%
%  The variables that are returned by this funtion are:
%        section... this is the section that is resampled.
%        x... This is the inline vector of the resampled section.  It will 
%                start at zero and continue to the end of the section.
%        z... This is the depth vector of the resampled section.  It will 
%                start at zero and continue to the end of the section.
%        dxz... This is the sampling rate that is used in both the depth 
%                and the inline directions.  ie. dz=dx=newdxz 
%
% By Heather J.E. Lloyd, August 2009
%
%NOTE: It is illegal for you to use this software for a purpose other
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

if nargin<3
    disp('At least oldsection, oldx and oldz must be supplied to use this function');
    return
end
if nargin==3
    dxz=10;
else
    dxz=newdxz;
end
szoldx=size(oldx);
if szoldx(1)>szoldx(2);
    nx=(oldx(1):dxz:oldx(end))';
    x=(0:dxz:(oldx(end)-oldx(1)))';
else
    x=0:dxz:(oldx(end)-oldx(1));
    nx=(oldx(1):dxz:oldx(end));
end

szoldz=size(oldz);
if szoldz(1)>szoldz(2);
    nz=(oldz(1):dxz:oldz(end))';
    z=(0:dxz:(oldz(end)-oldz(1)))';
else
    z=0:dxz:(oldz(end)-oldz(1));
    nz=(oldz(1):dxz:oldz(end));
end
section=interp2(oldx,oldz,oldsection,nx,nz);

end

