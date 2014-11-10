function [pressure_out]=ders_5pt(pressure,density,delx)
% DERS2_5PT ... compute the 5 point divergence (1/density * grad pressure)
%
% [pressure_out]=ders2_5pt(pressure,density,delx)
%
% DERS2_5PT computes the 5 point approximation of the spatial derivative 
% the two dimensional matrices 'input' and 'density.  The horizontal and
% vertical bin spacing of both matrices MUST be the same and equal.    
%
% pressure = input pressure matrix
% density = input density matrix
% delx = the horizontal/ vertical bin spacing in consistent units
%
% pressure_out = output pressure matrix
%
% by Carris Youzwishen, April 1999
% extended to full acoustic wave equation by Hugh Geiger, September 2003 
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

pres=zeros(nrows+2,ncolumns+2);
pres(2:nrows+1,2:ncolumns+1)=pressure;
pres(1,:)=0;
pres(end,:)=pres(end-1,:);
pres(:,1)=pres(:,2);
pres(:,end)=pres(:,end-1);
dens=zeros(nrows+2,ncolumns+2);
dens(2:nrows+1,2:ncolumns+1)=density;
dens(1,:)=dens(2,:);
dens(end,:)=dens(end-1,:);
dens(:,1)=dens(:,2);
dens(:,end)=dens(:,end-1);
factor=1/(2*delx^2);
clear pressure
clear density


% for first dimension

   pressure_out = (pres(3:nrows+2,2:ncolumns+1) - pres(2:nrows+1,2:ncolumns+1)).*...
                  (dens(3:nrows+2,2:ncolumns+1) + dens(2:nrows+1,2:ncolumns+1))./...
                   dens(3:nrows+2,2:ncolumns+1)*factor-...
                  (pres(2:nrows+1,2:ncolumns+1) - pres(1:nrows,2:ncolumns+1)).*...
                  (dens(2:nrows+1,2:ncolumns+1) + dens(1:nrows,2:ncolumns+1))./...
                   dens(1:nrows,2:ncolumns+1)*factor;

% for second dimension 

 pressure_out = (pres(2:nrows+1,3:ncolumns+2) - pres(2:nrows+1,2:ncolumns+1)).*...
                (dens(2:nrows+1,3:ncolumns+2) + dens(2:nrows+1,2:ncolumns+1))./...
                 dens(2:nrows+1,3:ncolumns+2)*factor-...
                (pres(2:nrows+1,2:ncolumns+1) - pres(2:nrows+1,1:ncolumns)).*...
                (dens(2:nrows+1,2:ncolumns+1) + dens(2:nrows+1,1:ncolumns))./...
                 dens(2:nrows+1,1:ncolumns)*factor+pressure_out;
