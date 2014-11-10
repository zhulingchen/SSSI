function [unow,uthen]=spongeABC(unow,uthen,nx,nz,nxsponge,nzsponge,coeff)
%
% spongeABC(unow,uthen,nx,nz,nxsponge,nzsponge,coeff)
%
% spongeABC: Sponge Absorbing boundary condtion
%
% 	unow:     the current wavefields
%
%	uthen:    the previous wavefields
%
%	nx:       model size of lateral direction
%
%	nz:       model size of vertical direction
%
%   nxsponge: sponge size of lateral direction
%
%   nzsponge: sponge size of vertical direction
%
%   coeff:    damping coefficient
%   Xiang Du, May   2007, $1.0
%             Sept. 2008, $1.1
%             Nov.  2008, $1.2
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

coeff2=coeff*coeff;

%left side
for t=1:nxsponge
  unow(:,t)=unow(:,t)*exp(-coeff2*(nxsponge-t)*(nxsponge-t));
  uthen(:,t)=uthen(:,t)*exp(-coeff2*(nxsponge-t)*(nxsponge-t));
end

%right side
for t=nx-nxsponge+1:nx
  unow(:,t)=unow(:,t)*exp(-coeff2*(t-nx+nxsponge-1)*(t-nx+nxsponge-1));
  uthen(:,t)=uthen(:,t)*exp(-coeff2*(t-nx+nxsponge-1)*(t-nx+nxsponge-1));
end

%bottom side
for t=nz-nzsponge:nz
  unow(t,:)=unow(t,:)*exp(-coeff2*(t-nz+nzsponge-1)*(t-nz+nzsponge-1));
  uthen(t,:)=uthen(t,:)*exp(-coeff2*(t-nz+nzsponge-1)*(t-nz+nzsponge-1));
end

%top side
for t=1:nzsponge
  unow(t,:)=unow(t,:)*exp(-coeff2*(nzsponge-t)*(nzsponge-t));
  uthen(t,:)=uthen(t,:)*exp(-coeff2*(nzsponge-t)*(nzsponge-t));
end
