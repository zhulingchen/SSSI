function a=alpine(m)
% a=alpine(m)
% 
if( nargin < 1) 
	m=size(get(gcf,'colormap'),1); 
end
ind =linspace(1,m,m);
g=((cos(2*pi*(ind-1)/(m-1))+1)/2)';
r=linspace(0.,1.,m)';
b = [linspace(0,1,m/2) ones(1,ceil(m/2))]';
a = [r.^2 g b.^2];
%a= brighten(a,-.5);
