function shotg=gaincc(shot,x,z,xshot,zshot)
% GAINCC ... gain a migrated shot record after CC imaging condition
%
% The method of gain correction is simply to compute the radial distance,
% r, from the shot to each image point and to multiply the migrated shot by
% r. For a constant velocity medium, this is the correct gain correction in
% 2D while for 3D it should be r^2. For variable velocity this is a kludge
% but seems to work reasonably well. The best way to gaincorrect a migrated
% shot is to use the deconvolution imaging condition.
%
% shot...input shot
% x ... x coordinates of each column of shot
% z ... z coordinates of each row of shot
% xshot ... xcoordinate of shot
% zshot ... z coordinate of shot
%
% G.F. Margrave, CREWES, 2010
%

%compute radial distance
[nz,nx]=size(shot);
x=x(:)';
z=z(:);
xx=ones(nz,1)*x;
zz=z*ones(1,nx);
r=sqrt((xx-xshot).^2+(zz-zshot).^2);

shotg=shot.*r;

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