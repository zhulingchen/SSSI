%  NORM_REFLECT.m
%
%  Calculates the normal incidence reflectivity matix given a velocity
%  model matrix.  May be modified to take in a density matrix as well.
%
%  refl=norm_reflect(vel,dens,clipvalue)
%
%  refl...........matrix of normal incidence reflectivities
%  vel............the velocity matrix
%  dens...........matrix of densities same size as vel(default = ones(size(vel))
%  clipvalue......the number of columns and rows to "clip" or remove from
%                 the edge of the matrix (for absorbing boundary
%                 conditions) (default=0)
%
%  Zoron Rodriguez, G.F. Margrave 2006
%
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.

% BEGIN TERMS OF USE LICENSE
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

function refl=norm_reflect(vel,dens,clipn)

if(nargin<3)
    clipn=0;
end
if (nargin<2);
    clipn=0;
    dens=ones(size(vel));
end

[nz,nx]=size(vel);

r=zeros(size(vel));

%loop over traces
for k=1:nx
    for s=2:nz
        I2=vel(s,k)*dens(s,k);
        I1=vel(s-1,k)*dens(s-1,k);
        r(s,k)= (I2-I1)/(I2+I1);
    end
end
   refl=zeros(size(r));
   refl(clipn+1:nz-clipn,clipn+1:nx-clipn)=r(clipn+1:nz-clipn,clipn+1:nx-clipn);
