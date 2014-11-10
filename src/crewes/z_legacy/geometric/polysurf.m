function c=polysurf(x,y,z,order,w)
% fit a 2-D polynomial surface to the data
%
% c=polysurf(x,y,z,order,w)
% c=polysurf(x,y,z,order)
%
% x = vector containing the x coordinates of the data
%			 
% y = vector containing the y coordinates of the data
%				  
% z = vector containing the z coordinates of the data
%				
% c= vector of 2-D polynomial coefficients
%
% w = vector of weights for the data points
%   ****** default = ones(size(z)) *******
%
% by G.F. Margrave, March 1993
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
%generate the 2-D Vandermonde matrix of order=order
numcols = sum( 1:order+1 );
V = zeros(length(z),numcols);
 
colcount = 1;
for n=order:-1:0
   for m=0:n
         V(:,colcount) = (x(:).^(n-m)).*(y(:).^m);
         colcount=colcount+1;
   end
end
% generate a matrix of weights and apply them to the Vandermonde matrix
if( nargin > 4 )
	[nr,nc]=size(V);
	[wr,wc]=size(w);
	if(wr==1) w=w';end
	z = z.*w;
	w=w*ones(1,nc);% each row gets a constant weight
	V=V.*w; % apply with array multiply
end
% now we have generated the matrix V necessary for the equation
% Vc=z where c is the unknown vector of 2-d polynomial coefficients
% The solution is:
% z=z.';
c = V\z;
 
 
