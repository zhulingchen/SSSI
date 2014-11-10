function plotfile(imagefile, nx, nz, dx, dz)

% plotfile: plot the migration file (binary format)
%
%   imagefile: migration file (binary format)
%
%	nx: the size of x direction
%
%   nz: the size of z direction
%
%   dx: the space interval of x direction
%
%   dz: the space interval of z direction
%
%   Xiang Du, Nov. 2008
%
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.

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

fid=fopen(imagefile,'r');
data=fread(fid,[nz nx],'float');
fclose(fid);

dmax = max(max(data));
dmin = min(min(data));
clim = [dmin dmax];

clim=clim*0.5;

x=[0:dx:(nx-1)*dx];
z=[0:dz:(nz-1)*dz];

imagesc(x,z, data, clim);
colormap(gray);
