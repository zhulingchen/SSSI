function yi=pcint(x,y,xi)
% PCINT: piecewise constant interpolation
%
% yi=pcint(x,y,xi)
%
% pcint does piecewise constant interpolation.
% That is, x and y are assumed to represent a function that is
% "blocky" or piecewise constant. This means that for any xi
% between x(k) and x(k+1), the function evaluates to y(k).
% Points with xi < x(1) evaluate to y(1) and for xi>x(length(x))
% evaluate to y(length(x))
% NOTE xi must be sorted in ascending order. (xi(k)<xi(k+1))
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

yi=nan*zeros(size(xi));
iiold=1;
nx=length(x);
%handle ends
j=find(xi<x(1));
if(~isempty(j))
	yi(j)=y(1)*ones(size(j));
end
j=find(xi>x(length(x)));
if(~isempty(j))
	yi(j)=y(length(x))*ones(size(j));
end
j=find(isnan(yi));
j=j(:)';
for k=j

		ii=find(x<=xi(k));

		yi(k)= y(ii(length(ii)));

end
