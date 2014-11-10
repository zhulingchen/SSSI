function [ix1,iz1,nxf,nzf,nxplot,initzp,nzplot].....
    = fdInitB3sizes(Dxz,nx,nz,mvXmax,mvZtop,mvZmax,wncMat)  %,iLbnd,shotX)
%Set up index arrays to be used in computations
%Allow for displacements on the left which simulate a 'mirrored surface'
%
%The input parameters are
%Dxz     .... FD sample rate in (metres)
%nx      .... Number of spatial samples in the X direction
%nz      .... Number of spatial samples in the Z direction
%mvXmax  .... X-length of movie frames to display or plot
%mvZtop  .... Z-level at top of each movie frame, often 0.
%mvZmax  .... Z-depth of movie frames to display or plot
%wncmat  .... A set of FD correction matrices, for particular Vp and Vs
%iLbnd   .... Boundary code left
%shotX   .... X (from the FD model left) of the initializing source
%
%The output parameters are
%ix1   ...... X co-ordinate of start of model in arrays
%iz1   ...... Z co-ordinate of start of model in arrays
%nxf     .... Number of X samples including border
%nzf     .... Number of Z samples including border
%nxplot  .... No. of X points to plot
%initzp  .... Initial Z point to plot
%nzplot  .... No. of Z points to plot
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

% nzm = nz-1;
% if iBrigid<=0
%     nzm = nz;
% end
% ixm = 1:nxf-1;
% ixm2 = 2:nxf-1;
% izm = 1:nzm;
% izm2 = 2:nzm;
% izb = nz+1;
% ivxi = ix1:ix1*2-2;     %Indicies to fill in symmetric displacements on left
% ivxo = ix1-1:-1:1;

%Allow for the extra displacement points needed for correction filtering

%disp(wncMat)
[nxC,nzC,nTerms] = size(wncMat);
%disp([nxC,nzC,nTerms])
nXextra = 0;
nZextra = 0;
if nTerms == 6      %Indicator of legitimate correction operator
    nXextra = (nxC-1)/2;    %Not necessary for interior source
    nZextra = (nzC-1)/2;
end
%     iExtrai = 1:nExtra+1;
%     iExtrao = -1:-1:-(nExtra+1);
%     iExtra = 1:nExtra;
%     zFill = zeros(1,nExtra);
%disp(nExtra)
% ix1 = nXextra + 1;
% iz1 = nZextra + 1;
% nxf = nx + nXextra*2;
% nzf = nz + nZextra*2;
ix1 = nXextra + 2;          %One extra column for the FD
iz1 = nZextra + 2;
% nxf = nx + nXextra*2+2;
% nzf = nz + nZextra*2+2;
nxf = nx + nXextra*2 + 2;
nzf = nz + nZextra*2 + 2;
nxplot = round(mvXmax/Dxz)+ix1; 
initzp = round(mvZtop/Dxz)+iz1; 
nzplot = round(mvZmax/Dxz)+iz1;
%Will be used as follows
    %Uz(ix1+iExtrao) = Uz(ix1+iExtrai);     %Symmetric edge
    %Ux(ix1+iExtrao) = Ux(ix1+iExtrai-1);   %Symmetric edge
    %Uz(nx+iExtra+1) = zFill;               %Rigid edge
    %Ux(nx+iExtra+1) = zFill;               %Rigid edge

