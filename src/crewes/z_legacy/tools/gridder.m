function [xi,yi,zi] = gridder(x,y,z,xi,yi, bndyval, colbndyfac, rowbndyfac,...
			bndydist)
% GRIDDER Data gridding.
%
% [xi,yi,zi] = gridder(x,y,z,xi,yi, bndyval, colbndyfac, rowbndyfac, bndydist)
% [xi,yi,zi] = gridder(x,y,z,xi,yi, bndyval, colbndyfac, rowbndyfac)
% [xi,yi,zi] = gridder(x,y,z,xi,yi, bndyval, colbndyfac)
% [xi,yi,zi] = gridder(x,y,z,xi,yi, bndyval)
% [xi,yi,zi] = gridder(x,y,z,xi,yi)
%
%	x =  vector of x coordinates of the random points to be gridded
%	y =  vector of y coordinates of the random points to be gridded
%	z =  vector of z values of the random points to be gridded
%	**** requirement: x,y,z must all be the same length ****
%       xi = a vector or matrix of x coordinates describing the output
%		grid. If a vector, xi describes the x coordinates of 
%               each column of the grid.
%       yi = a vector or matrix of y coordinates describing the output
%		grid. If a vector, yi describes the y coordinates of 
%               each row of the grid.
%	bndyval = scalar giving the boundary value
%	**** default is no boundary value constraints *****
%	colbndyfac = scalar specifing the number of boundary values to 
%		be seeded along the top and bottom of the columns. One 
%		boundary value will be placed for every colbndyfac'th
%		column.
%       **** default = 10 ****
%	rowbndyfac = scalar specifying the number of row boundary values.
%		One boundary value will be place for every rowbndyfac'th
%		row.
%       **** default = 10 ***
%	bndydist = 'rand' ... distribute boundary values randomly
%		   'reg'  ... distribute boundary values regularly
%       *** default = 'reg' ***
%
%	GRIDDER is a modification of Matlab's GRIDDATA algorithm with
%	two major changes:
%		1) Matlab's routine fails if any of the input points have
%		the same (x,y) coordinates. This was resolved by solving
%		for the data weights using a pseudo matrix inverse instead
%		of the hard matrix division that GRIDDATA uses. The result 
%		is a much more stable algorithm at the cost of slightly 
%		greater run times.
%		2) GRIDDER allows the specification of 'boundary values' 
%		for the grid to control edge effects. This is done by 
%		augmenting the input points (x,y,z) with a set of values
%		(xb,yb,zb) spread around the outside of the grid (xi,yi).
%		This is very effective in controling the edges but causes
%		an increase in run time which is quadratic with the number of
%		boundary points.
%
% Changes by: G.F. Margrave, October 1993
%
%	ZI = GRIDDATA(X,Y,Z,XI,YI) returns matrix ZI containing elements
%	corresponding to the elements of matrices XI and YI and determined by
%	interpolation within the 2-D function described by the (usually)
%	nonuniformly-spaced vectors (X,Y,Z).
%
%	XI can be a row vector, in which case it specifies a matrix with
%	constant columns. Similarly, YI can be a column vector and it 
%	specifies a matrix with constant rows. 
%
%	[XI,YI,ZI] = GRIDDATA(X,Y,Z,XI,YI) returns the XI and YI formed
%	this way, which are the same as the matrices returned by MESHGRID.
%
%	GRIDDATA uses an inverse distance method.
%
%	See also INTERP2, INTERP1.
%	Copyright (c) 1984-93 by The MathWorks, Inc.
%       Reference:  David T. Sandwell, Biharmonic spline
%       interpolation of GEOS-3 and SEASAT altimeter
%       data, Geophysical Research Letters, 2, 139-142,
%       1987.  Describes interpolation using value or
%       gradient of value in any dimension.
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
error(nargchk(5,9,nargin)) % Ignore method for now.
[msg,x,y,z,xi,yi] = xyzchk(x,y,z,xi,yi);
if length(msg)>0, error(msg); end
% generate the boundary values
   if( nargin > 5 )
	dely = yi(2,1)-yi(1,1);
	delx = xi(1,2)-xi(1,1);
	[m,n]=size(xi);
	if( nargin < 7) colbndyfac = 10; end
	if( nargin < 8) rowbndyfac = 10; end
	if( nargin < 9) bndydist = 'reg';end
	if( strcmp( bndydist,'reg' ) )
		% top values
		ntop = n/colbndyfac;
		ytop = (yi(1,1)-dely)*ones(1,ntop);
		xtop = linspace(xi(1,1),xi(1,n),ntop);
		ztop = bndyval*ones(1,ntop);
		% bottom values
		nbottom = n/colbndyfac;
		ybottom = (yi(m,1)+dely)*ones(1,nbottom);
		xbottom = xtop;
		zbottom = ztop;
		% left values
		nleft = m/rowbndyfac;
		xleft = (xi(1,1)-delx)*ones(1,nleft);
		yleft = linspace(yi(1,1),yi(m,1),nleft);
		zleft = bndyval*ones(1,nleft);
		% right values
		xright = (xi(1,n)+delx)*ones(1,nleft);
		yright = linspace(yi(1,1),yi(m,1),nleft);
		zright = zleft;
	elseif( strcmp( bndydist,'rand') )
		% top values
		ntop = n/colbndyfac;
		y1 = yi(1,1); y2 = y1-2*dely;
		ytop = rand(1,ntop)*(y2-y1)+y1;
		x1=xi(1,1);x2=xi(1,n);
		xtop = rand(1,ntop)*(x2-x1)+x1;
		ztop = bndyval*ones(1,ntop);
		%bottom values
		y1 = yi(m,1); y2 = y1+2*dely;
		ybottom = rand(1,ntop)*(y2-y1)+y1;
		x1=xi(1,1);x2=xi(1,n);
		xbottom = rand(1,ntop)*(x2-x1)+x1;
		zbottom = bndyval*ones(1,ntop);
		%left values
		nleft = m/rowbndyfac;
		y1 = yi(1,1); y2 = yi(m,1);
		yleft = rand(1,nleft)*(y2-y1)+y1;
		x1=xi(1,1);x2=x1-2*delx;
		xleft = rand(1,nleft)*(x2-x1)+x1;
		zleft = bndyval*ones(1,nleft);
		%right values
		y1 = yi(1,1); y2 = yi(m,1);
		yright = rand(1,nleft)*(y2-y1)+y1;
		x1=xi(1,n);x2=x1+2*delx;
		xright = rand(1,nleft)*(x2-x1)+x1;
		zright = bndyval*ones(1,nleft);
	else
		error('illegal value for ''bndydist'' ');
	end
	% attach the values on the end of the input points;
	x = [x(:); xtop(:); xbottom(:); xleft(:); xright(:)];
	y = [y(:); ytop(:); ybottom(:); yleft(:); yright(:)];
	z = [z(:); ztop(:); zbottom(:); zleft(:); zright(:)];
   end
		
xy = x(:) + y(:)*sqrt(-1);
		
% Determine weights for interpolation
d = xy * ones(1,length(xy));
d = abs(d - d.');
mask = find(d == 0);
d(mask) = ones(length(mask),1);
g = (d.^2) .* (log(d)-1);   % Green's function.
g(mask) = zeros(length(mask),1); % Value of Green's function at zero
weights = pinv(g)*z(:);
[m,n] = size(xi);
zi = zeros(m,n);
jay = sqrt(-1);
xy = xy.';
% Evaluate at requested points (xi,yi).  Loop to save memory.
for i=1:m
  for j=1:n
    d = abs(xi(i,j)+jay*yi(i,j) - xy);
    mask = find(d == 0);
    if length(mask)>0, d(mask) = ones(length(mask),1); end
    g = (d.^2) .* (log(d)-1);   % Green's function.
    % Value of Green's function at zero
    if length(mask)>0, g(mask) = zeros(length(mask),1); end
    zi(i,j) = g * weights;
  end
end
if nargout<=1,
  xi = zi;
end
