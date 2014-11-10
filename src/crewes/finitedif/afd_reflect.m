function r=afd_reflect(velocity,clipn)
% AFD_REFLECT ... calculate the reflectivity from a velocity model
%
% r=afd_reflect(velocity,clipn)
% 
% AFD_REFLECT calculates the normal indicence reflectivity of an input
% velocity model for use with AFD_EXPLOAD.  
%
% velocity = velocity matrix
% clipn = the number of "bin layers" with which to remove from the
%       edge of velocity matrix in order to prevent artifacts from
%       the absorbing boundary conditions
%     = the suggested number is 5, but zero may also be chosen
% r = refelctivity calculated as .5*abs(gradient(log(velocity)))
%
% Carrie Youzwishen, April 1999
% completely rewritten by G.F. Margrave August 2000
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

% recalculate velocity for exploding reflector model (1 way travel time)

[nz,nx]=size(velocity);

[rx,rz]=gradient(log(velocity));

tmp=.5*sign(rz).*sqrt(rx.^2+rz.^2);
clear rx rz

r=zeros(size(tmp));
r(clipn+1:nz-clipn,clipn+1:nx-clipn)=tmp(clipn+1:nz-clipn,clipn+1:nx-clipn);





%create reflectivity matrix
%vnew=zeros(nz+2,nx);
%r=zeros(nz,nx);

%vnew(2:nz+1,:) = 1/2*log(velocity);
%vnew(1,:)=vnew(2,:);
%vnew(nz+2,:)=vnew(nz+1,:);

%temp(:,:)=(vnew(3:nz+2,:) - vnew(1:nz,:))/2;
%r(clipn+1:nz-clipn,clipn+1:nx-clipn)=temp(clipn+1:nz-clipn,clipn+1:nx-clipn);
