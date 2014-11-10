function [pressure_out]=spatial_derivs_order2(pressure,logdensity,delx)
% spatial_derivs_order2 ... second order spatial derivatives for acoustic propagator
%
% [pressure_out]=spatial_derivs_order2(pressure,logdensity,delx)
%
% spatial_derivs_order2 computes the second order approximation of the spatial
% derivative term of the acoustic wave equation.  This is
% Laplacian(pressure)-grad(logdensity)_dot_grad(pressure).  The horizontal
% and vertical bin spacing of both matrices MUST be the same and equal.
%
% pressure = input pressure matrix
% logdensity = input logdensity matrix
% delx = the horizontal/ vertical bin spacing in consistent units
%
% pressure_out = output pressure matrix
%
% by Carris Youzwishen, April 1999
% extended to full acoustic wave equation by Gary Margrave, Feb 2014
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

[nrows,ncolumns]=size(pressure);

pres=zeros(nrows+2,ncolumns+2);%new matrix with an extra row and column all around
pres(2:nrows+1,2:ncolumns+1)=pressure;%throw pressure into the middle
% pres(1,:)=pres(2,:);%fill in extra row on top
%NOTE: Flooding the top extra row with live data values causes a very
%strong surface wave. Better to leave it zero. I think this may have
%something to do with boundary conditions.
pres(end,:)=pres(end-1,:);%fill extra row on bottom
pres(:,1)=pres(:,2);%fill extra column on left
pres(:,end)=pres(:,end-1);%fill extra column on right
dens=zeros(nrows+2,ncolumns+2);%new density matrix
dens(2:nrows+1,2:ncolumns+1)=logdensity;%fill
dens(1,:)=dens(2,:);%extra row on top
dens(end,:)=dens(end-1,:);%extra row on bottom
dens(:,1)=dens(:,2);%extra column on left
dens(:,end)=dens(:,end-1);%extra column on right
factor=1/(delx^2);%precompute factor
clear pressure
clear logdensity


%laplacian of pressure
% for first dimension

   pressure_out = (pres(3:nrows+2,2:ncolumns+1) - 2*pres(2:nrows+1,2:ncolumns+1) +...
                  pres(1:nrows,2:ncolumns+1))*factor;

% for second dimension 

   pressure_out = pressure_out + (pres(2:nrows+1,3:ncolumns+2) - 2*pres(2:nrows+1,2:ncolumns+1) +...
                 pres(2:nrows+1,1:ncolumns))*factor;
             
%   output = (input2(3:nrows+2,2:ncolumns+1) - 2.*input2(2:nrows+1,2:ncolumns+1) +...
%                   input2(1:nrows,2:ncolumns+1))*factor;
% 
% % for second dimension 
% 
%  output = (input2(2:nrows+1,3:ncolumns+2) - 2.*input2(2:nrows+1,2:ncolumns+1) +...
%                  input2(2:nrows+1,1:ncolumns))*factor+ output;             
             
%grad(logdensity)_dot_grad(pressure)
   factor=.25*factor;
%first dimension
   pressure_out=pressure_out - ...
       ((dens(3:nrows+2,2:ncolumns+1)-dens(1:nrows, 2:ncolumns+1)).*...
       (pres(3:nrows+2,2:ncolumns+1)-pres(1:nrows, 2:ncolumns+1)))*factor;
%second dimension
   pressure_out=pressure_out - ...
       ((dens(2:nrows+1,3:ncolumns+2)-dens(2:nrows+1, 1:ncolumns)).*...
       (pres(2:nrows+1,3:ncolumns+2)-pres(2:nrows+1, 1:ncolumns)))*factor;

