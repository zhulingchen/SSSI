function [xvect,yvect,zvect] = gen3Arr(nx,ny,nz,xzero,yzero,zzero)
%function [xvect,yvect,zvect] = gen3Arr(nx,ny,nz,xzero,yzero,zzero)
%Generate x, y, and z vectors that will form an (nx,ny,nz) array, starting at
%   xzero,yzero,zzero and incrementing by 1
%The input parameters are
%nx      .... No. of points in the X-direction
%ny      .... No. of points in the Y-direction
%nz      .... No. of points in the Z-direction
%xzero   .... Starting point of the X-aray
%yzero   .... Starting point of the Y-aray
%zzero   .... Starting point of the Z-aray
%The output parameters are
%xvect   .... For any no. n in range
%yvect   .... xvect(n), yvect(n), zvect(n)
%zvect   .... is a unique point
%
% P.M. Manning, Dec 2011
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
    %xzero, zzero       Increments are all 1
xvect = zeros(1,nx*ny*nz);
yvect = xvect;
zvect = xvect;
ind = 1;
z = zzero;
for iz = 1:nz
    y = yzero;
    for iy = 1:ny
        x = xzero;
        for ix = 1:nx
            xvect(ind) = x;
            yvect(ind) = y;
            zvect(ind) = z;
            ind = ind+1;
            x = x+1;
        end
        y = y+1;
    end
    z = z+1;
end