function s=seisclrs_old(m,pct)
% SEISCLRS: creates the default color map for plotimage.
%
% s=seisclrs(m,pct)
%
% seisclrs creates a gray-level color map designed to display
% seismic data using plotimage.
% m ... number of gray levels
%  ***** default is determined by the colormap size in gcf *****
% pct ... percentage of the colormap that transitions between 
%	black and white. 100 gives a linear ramp while 1 gives an
%	abrupt transition.
%  ***** default is 50 *******
% s ... matix of size [m,3] giving the gray levels. Each column is
%	is identical.
%
if(nargin < 2) pct=50; end
if( nargin < 1) m=size(get(gcf,'colormap'),1); end
if(pct>100 | pct <1) disp('gray_percentage set poorly in seisclrs'); end
glen = round(.5*pct*m/100);
blen=floor(m/2)-glen;
wlen=m-blen-2*glen;
%blen=m/2-glen;
%wlen=m-blen-2*glen;
grad=linspace(0,1,2*glen);
g=[zeros(blen,1);grad(1:2*glen)';ones(wlen,1)];
g=flipud(g);
s = [g g g];
global NOBRIGHTEN
if(isempty(NOBRIGHTEN)) nobrighten=0; else nobrighten=NOBRIGHTEN; end
if(~nobrighten)
	s= brighten(s,.5);
end
