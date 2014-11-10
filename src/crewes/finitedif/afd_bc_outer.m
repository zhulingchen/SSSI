function [snapshot]=afd_bc_outer(delx,delt,velocity,snap1,snap2,snapshot,boundary)
% AFD_BC_OUTER ... apply absorbing boundary condition to outer boundary
%
% [snapshot]=afd_bc_outer(delx,delt,velocity,snap1,snap2,snapshot,boundary)
% AFD_BC_OUTER applies absorbing boundary conditions to the outer layer of rows
% and columns in a situation where there are 1 or 2 layers of rows and columns on
% the boundary.  The output matrix has absorbing boundary conditions along its
% outer edge. In order to be stable, vmax*delt/delx must be < 0.7.  
%
% delx = the grid spacing for horizontal and vertical bins in consistent units
% delt = the time interval in seconds
% velocity = the input velocity matrix in consisnent units
%          = has a size of floor(zmax/delx)+1 by floor(xmax/delx)+1
% snap1 = the wavefield at time=t-2*delt (same size as velocity matrix)
%        = will be based on the source array desired i.e. the position
%          of the sources will be one, and the rest of the positions
%          will be zero
% snap2 = the wavefield at time = t-delt(same size as velocity matrix)
% snapshot = the wavefield at time = t (same size as velocity matrix)
% boundary = indicate whether all sides of the matrix are absorbing
%          = 0 indicates that no absorbing boundaries are desired
%          = 1 indicates all four sides are absorbing
%          = 2 choses three sides to be absorbing, and the top one not to be
%             this enables sources to be put on the surface
%
% snapshot = the wavefield at time = t with absorbing boundary conditions
%
% by Carrie Youzwishen, Febraury 1999 
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

[nz,nx]=size(snap1);

if boundary == 1

  % apply absorbing boundary for top 

  snapshot(1,3:nx-2) = (2.*velocity(1,3:nx-2).*delx*delt.^2)./(delx+velocity(1,3:nx-2).*delt).*...
         ((snapshot(2,3:nx-2)./(2*delt*delx) - snap1(2,3:nx-2)./(2*delx*delt) + ....
         snap1(1,3:nx-2)./(2*delt*delx)) + ...
    (-1./(2.*delt.^2.*velocity(1,3:nx-2))).*...
         (-2.*snap2(1,3:nx-2) + snap1(1,3:nx-2 ) -2.*snap2(2,3:nx-2)+ ...
         snap1(2,3:nx-2) + snapshot(2,3:nx-2)) + ...
    velocity(1,3:nx-2)./(4*delx.^2).*...
         (snapshot(2,4:nx-1) + snap1(1,4:nx-1) + snapshot(2,2:nx-3) - ...
          2.* snapshot(2,3:nx-2) - 2.*snap1(1,3:nx-2) + snap1(1,2:nx-3)));
end

% apply absorbing boundary for bottom

snapshot(nz,3:nx-2) = -2.*delx.*delt.^2.*velocity(nz,3:nx-2)./(delx+velocity(nz,3:nx-2).*delt).*...
       ((-snapshot(nz-1,3:nx-2)./(2*delt*delx) - snap1(nz,3:nx-2)./(2*delt*delx) +...
       snap1(nz-1,3:nx-2)./(2*delt*delx)) +...
  1./(2*delt.^2.*velocity(nz,3:nx-2)).*...
       (-2.*snap2(nz,3:nx-2) + snap1(nz,3:nx-2) + snapshot(nz-1,3:nx-2) - ...
       2.*snap2(nz-1,3:nx-2) + snap1(nz-1,3:nx-2)) + ...
  (-velocity(nz,3:nx-2)./(4*delx.^2)).* ...
       (snapshot(nz-1,4:nx-1) - 2.*snapshot(nz-1,3:nx-2) + snapshot(nz-1,2:nx-3) + ...
       snap1(nz,4:nx-1) - 2.*snap1(nz,3:nx-2) + snap1(nz,2:nx-3)));

% apply absorbing boundary for right hand side

snapshot(3:nz-2,nx) =  -2.*delx.*delt.^2.*velocity(3:nz-2,nx)./(delx+velocity(3:nz-2,nx).*delt).*...
       ((-snapshot(3:nz-2,nx-1)./(2*delt*delx) - snap1(3:nz-2,nx)./(2*delt*delx) + ...
       snap1(3:nz-2,nx-1)./(2*delt*delx)) + ...
  1./(2*delt.^2.*velocity(3:nz-2,nx)).*...
       (-2.*snap2(3:nz-2,nx) + snap1(3:nz-2,nx) + snapshot(3:nz-2,nx-1) - ...
       2.*snap2(3:nz-2,nx-1) + snap1(3:nz-2,nx-1)) + ...
  (-velocity(3:nz-2,nx)./(4*delx.^2)).* ...
       (snapshot(4:nz-1,nx-1) - 2.*snapshot(3:nz-2,nx-1) + snapshot(2:nz-3,nx-1) + ...
       snap1(4:nz-1,nx) - 2.*snap1(3:nz-2,nx) + snap1(2:nz-3,nx)));


% apply absorbing boundary for left hand side

snapshot(3:nz-2,1) =(2.*velocity(3:nz-2,1).*delx*delt.^2)./(delx+velocity(3:nz-2,1).*delt).*...
       ((snapshot(3:nz-2,2)./(2*delt*delx) - snap1(3:nz-2,2)./(2*delt*delx) + ...
       snap1(3:nz-2,1)./(2*delt*delx)) + ...
  (-1./(2*delt.^2.*velocity(3:nz-2,1))).* ...
       (-2.*snap2(3:nz-2,1) + snap1(3:nz-2,1) + snapshot(3:nz-2,2) - ...
       2.*snap2(3:nz-2,2) + snap1(3:nz-2,2)) + ...
  (velocity(3:nz-2,1)./(4*delx.^2)).* ...
       (snapshot(4:nz-1,2) - 2.*snapshot(3:nz-2,2) + snapshot(2:nz-3,2) +...
       snap1(4:nz-1,1) - 2.*snap1(3:nz-2,1) + snap1(2:nz-3,1)));


%%%%%% Absorbing boundary conditions for corners

% for lower right hand corner

snapshot(nz-1,nx) = velocity(nz-1,nx).*delt.*delx./(2.*velocity(nz-1,nx).*delt + (2^1/2)*delx).*...
       (snapshot(nz-2,nx)./delx + snapshot(nz-1,nx-1)./delx + ...
       (2^1/2)/(velocity(nz-1,nx).*delt)*snap2(nz-1,nx));

snapshot(nz,nx-1) = velocity(nz,nx-1).*delt.*delx./(2.*velocity(nz,nx-1).*delt + (2^1/2)*delx).*...
       (snapshot(nz-1,nx-1)./delx + snapshot(nz,nx-2)./delx + ...
       (2^1/2)/(velocity(nz,nx-1).*delt)*snap2(nz,nx-1));

snapshot(nz,nx) = velocity(nz,nx).*delt.*delx./(2.*velocity(nz,nx).*delt + (2^1/2)*delx).*...
       (snapshot(nz-1,nx)./delx + snapshot(nz,nx-1)./delx + ...
       (2^1/2)/(velocity(nz,nx).*delt)*snap2(nz,nx));

% for lower left hand corner

snapshot(nz-1,1) = velocity(nz-1,1).*delt.*delx./(2.*velocity(nz-1,1).*delt + (2^1/2)*delx).*...
       (snapshot(nz-2,1)/delx + snapshot(nz-1,2)/delx +...
       (2^1/2)/(velocity(nz-1,1)*delt)*snap2(nz-1,1));

snapshot(nz,2) = velocity(nz,2).*delt.*delx./(2.*velocity(nz,2).*delt + (2^1/2)*delx).*...
       (snapshot(nz-1,2)/delx + snapshot(nz,3)/delx +...
       (2^1/2)/(velocity(nz,2)*delt)*snap2(nz,2));

snapshot(nz,1) = velocity(nz,1).*delt.*delx./(2.*velocity(nz,1).*delt + (2^1/2)*delx).*...
       (snapshot(nz-1,1)/delx + snapshot(nz,2)/delx +...
       (2^1/2)/(velocity(nz,1)*delt)*snap2(nz,1));

if boundary == 1

  % for upper right hand corner

  snapshot(2,nx) = velocity(2,nx).*delt.*delx./(2.*velocity(2,nx).*delt + (2^1/2)*delx).*...
         (snapshot(3,nx)/delx + snapshot(2,nx-1)/delx +...
         (2^1/2)/(velocity(2,nx)*delt)*snap2(2,nx));

  snapshot(1,nx-1) = velocity(1,nx-1).*delt.*delx./(2.*velocity(1,nx-1).*delt + (2^1/2)*delx).*...
         (snapshot(2,nx-1)/delx + snapshot(1,nx-2)/delx +...
         (2^1/2)/(velocity(1,nx-1)*delt)*snap2(1,nx-1));

  snapshot(1,nx) = velocity(1,nx).*delt.*delx./(2.*velocity(1,nx).*delt + (2^1/2)*delx).*...
         (snapshot(2,nx)/delx + snapshot(1,nx-1)/delx +...
         (2^1/2)/(velocity(1,nx)*delt)*snap2(1,nx));

  % for upper left hand corner

  snapshot(2,1) = velocity(2,1).*delt.*delx./(2.*velocity(2,1).*delt + (2^1/2)*delx).*...
         (snapshot(3,1)/delx + snapshot(2,2)/delx + ...
         (2^1/2)/(velocity(2,1)*delt)*snap2(2,1));
 
  snapshot(1,2) = velocity(1,2).*delt.*delx./(2.*velocity(1,2).*delt + (2^1/2)*delx).*...
         (snapshot(2,2)/delx + snapshot(1,3)/delx + ...
         (2^1/2)/(velocity(1,2)*delt)*snap2(1,2));

  snapshot(1,1) = velocity(1,1).*delt.*delx./(2.*velocity(1,1).*delt + (2^1/2)*delx).*...
         (snapshot(2,1)/delx + snapshot(1,2)/delx + ...
         (2^1/2)/(velocity(1,1)*delt)*snap2(1,1));

end



