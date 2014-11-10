function s=seis(m,pct)
if(nargin < 2) pct=20; end
if( nargin < 1) m=size(get(gcf,'colormap'),1); end
if(pct>100 | pct <1) disp('gray_percentage set poorly in seisclrs'); end
glen = round(.5*pct*m/100);
blen=floor(m/2)-glen;
wlen=m-blen-2*glen;
%blen=m/2-glen;
%wlen=m-blen-2*glen;
grad=linspace(0,1,2*glen);
g=[zeros(blen,1);grad';ones(wlen,1)];
g=flipud(g);
s = [g g g];
s= brighten(s,.5);
