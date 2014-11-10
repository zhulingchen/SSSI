function h=drawvint(z,v,kol,ls)
% DRAWVINT: draw interval velocity as a piecewise constant function
%
% h=drawvint(z,v,kol,ls)
%
% draw an interval velocity function
%
% v= interval velocity vector
% z= depth vector
% kol = color
% default 'b'
% ls = linestyle to use
% default '-'
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


hax=gca;

if(nargin<4)
	ls='-';
end

if(nargin<3)
	kol='b';
end

dz=mean(diff(z));

z=z(:);
v=v(:);
z2=zeros(2*length(z),1);
z2(1:2:length(z2))=z;
z2(2:2:length(z2)-1)=z(2:length(z));
z2(length(z2))=z(length(z))+dz;

v2=zeros(2*length(z),1);
v2(1:2:length(v2))=v;
v2(2:2:length(v2))=v;

h=line(v2,z2,'color',kol,'linestyle',ls);

