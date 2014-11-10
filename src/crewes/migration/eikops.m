function [trM]=eikops(locs,z,velmod,xv,zv)
% EIKOPS: Eikonal Time Calculation - For PS (Converted) Wave migration.
% 
% [trM]=eikops(locs,x,z,velmod,xv,zv);
%
% locs... subindex in xv of the locations to be calculated 
% z .... elevation for each xv location (receiver) - a vector 
% velmod ... velocity model. Can be P wave or S wave. This is a matrix of 
%           velocities as a function of lateral position and depth.
% xv ... space coordinate vector for velmod
%           Requirement: length(xv)=size(velmod,2);
% zv ... space coordinate vector for velmod
%           Requirement: length(zv)=size(velmod,1);
%
% OUTPUT arguments:
%    trM ...all locations traveltimes: 3-D matrix with (P or S wave) traveltimes 
%            for the geological model.
%            Dimensions: zv*xv*nr, where nr = number of locations.


% EIKOPS calculates wave propagation time tables for Converted (PS) wave single shot
%   in Kirchhoff prestack depth migration. The algorithm used is the implementation of 
%   the Eikonal method, as presented by Sethian and Popovici, 1999: "3-D traveltime 
%   computation using the fast marching method" (Geophysics, 64, 516-523),
%   which was coded by Chad Hogan at CREWES.
%   The input velocity 'velmod' defines if it is an S- or P-wave
%
% Adapted to be used with a PS-wave Kirchhoff PSDM from topography by Saul
% Guevara
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

%tstart=clock;
[nvsamp,nvtr]=size(velmod);

% check the validity input arguments

dx=xv(2)-xv(1);
dz=dx;

%  ---- test velocity info ----

if(length(zv)~=nvsamp)
    error('Depth vector for velocity model is incorrect')   
end
if(length(xv)~=nvtr)
    error('Space coordinate vector for velocity model is incorrect')
end


% Time control        
clock1=cputime;
ntimes=0;
ievery=5;%print a progress message every this many traces


disp(['Starting time tables (3D matrix): '])
        
% Location Time: Wave propagation from all the locations (assuming as if they
% were sources)

        nr=length(locs);
        nz=length(zv);nx=length(xv);
        trM=zeros(nz,nx,nr);
        xr=xv(locs);
        zr=z(locs);
        
     % % Place for a code to edit zero velocities (not permitted for the 
     % %       Eikonal code)
     for i=1:nz,for j=1:nx,if velmod(i,j)==0,velmod(i,j)=0.1;end;end;end
    keyboard  
    for ikr=1:nr,
        xloc=round(xr(ikr)/dx)+1;
        zloc=round(zr(ikr)/dz)+1;        
        disp([' Calculating table for location no.',int2str(ikr),' of ' int2str(nr)])
        %disp(['Receiver: ',int2str(xloc), ' of ',int2str(length(x)) ])
        trtemp=eikonal2D(velmod,dx,xloc,zloc);

        if(rem(ikr,ievery)==0)
	    disp([' Completed table no. ' ,int2str(ikr) ,' of ' int2str(nr) ]);
        end
        
     % % Place for a code to edit the very long times due to the low velocity
     % %        (It is assumed that anything above 10 s is a long time)
    for i=1:nz,
        for j=1:nx,
            if trtemp(i,j)>10,trtemp(i,j)=0;
            end;
        end;
    end;

      
        trM(:,:,ikr)=trtemp ;
    end
    
totaltime=cputime-clock1;
disp(['Total time required ' num2str(totaltime)])



