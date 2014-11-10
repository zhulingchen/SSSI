function [snapshot,z,x]=afd_snap(delx,delt,velocity,snap1,snap2,laplacian,boundary)
% AFD_SNAP ... take one finite difference time step
%
% [snapshot,z,x]=afd_snap(delx,delt,velocity,snap1,snap2,laplacian,boundary)
%
% AFD_SNAP propogates a wavefield forward in depth by 
% one time step.  Two input matrices of the wavefield, one at 
% time=0-delt and one at time=0, are used in a finite 
% difference algorithm to propogate the wavefield.  The 
% finite difference algorithm can be calculated with a 
% five or nine point approximation to the Laplacian operator.  
% The five point approximation is faster, but the nine 
% point results in a broader bandwidth.The snapshot of this 
% propagated wavefield is returned. Note that the velocity 
% and grid spacing must fulfill the equation
% max(velocity)*delt/delx > 0.7 for the model to be stable.  
% This condition usually results in snap1 and snap2 
% being identical.    
%
% delx = the horizontal AND vertical bin spacing in consistent units
% delt = time interval in seconds
% velocity = the input velocity matrix in consisnent units
%          = has a size of floor(zmax/delx)+1 by floor(xmax/delx)+1
% snap1 = the wavefield at time=0 - delt (same size as velocity matrix)
%        = will be based on the source array desired i.e. the position
%          of the sources will be one, and the rest of the positions
%          will be zero
% snap2 = the wavefield at time = 0 (same size as velocity matrix)
% laplacian = an option between two approximations to the laplacian operator
%           = 1 is a 5 point approximation
%           = 2 is a nine point approximation
% boundary = indicate whether all sides of the matrix are absorbing
%          = 0 indicates that no absorbing boundaries are desired
%          = 1 indicates all four sides are absorbing
%          = 2 choses three sides to be absorbing, and the top one not to be
%             this enables sources to be put on the surface
%
% snapshot = the wavefield propagated forward one time interval
%            where the time interval = delt
% 
% by Carrie Youzwishen, February 1999
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
if(prod(size(snap1))~=prod(size(snap2)))
	error('snap1 and snap2 must be the same size');
end
xmax=(nx-1)*delx;
zmax=(nz-1)*delx;

x=0:delx:xmax;
z=(0:delx:zmax)';

if laplacian == 2

  snapshot=(velocity.^2.*delt^2).*del2_9pt(snap2,delx) + 2*snap2 - snap1;

   %prepare for absorbing bc's by zeroing outer 2 rows and columns
   if boundary == 1
   	  snapshot(1:2,:)=zeros(2,nx);
   	  snapshot(nz-1:nz,:)=zeros(2,nx);
   	  snapshot(:,1:2)=zeros(nz,2);
   	  snapshot(:,nx-1:nx)=zeros(nz,2);
   else
      snapshot(nz-1:nz,:)=zeros(2,nx);
   	  snapshot(:,1:2)=zeros(nz,2);
   	  snapshot(:,nx-1:nx)=zeros(nz,2);
   end

   if(boundary)
      [snapshot]=afd_bc_inner(delx,delt,velocity,snap1,snap2,snapshot,boundary);

      [snapshot]=afd_bc_outer(delx,delt,velocity,snap1,snap2,snapshot,boundary);
   end

else

  snapshot=velocity.^2.*delt^2.*del2_5pt(snap2,delx) + 2*snap2 - snap1;
	
   %prepare for absorbing bc's by zeroing outer 1 row and column
   if boundary == 1
     snapshot(1,:)=zeros(1,nx);
   	 snapshot(nz,:)=zeros(1,nx);
   	 snapshot(:,1)=zeros(nz,1);
   	 snapshot(:,nx)=zeros(nz,1);
   else
     snapshot(nz,:)=zeros(1,nx);
   	 snapshot(:,1)=zeros(nz,1);
   	 snapshot(:,nx)=zeros(nz,1);
   end

   if(boundary)
      [snapshot]=afd_bc_outer(delx,delt,velocity,snap1,snap2,snapshot,boundary);
   end  

end


