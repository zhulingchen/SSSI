function plotwigl2(amat,t,x,tstart,tend)

%function plotwigl(amat,t,x)
%plot wiggle trace from standard format

nx = length(x);
s = size(amat);
vmax = max(max(abs(amat)));
xmax = max(x);
scale = xmax/(nx*vmax);

amat = amat.*scale*3;

for k=1:nx
   amat(:,k) = amat(:,k) + x(k);
end

if(nargin==5)
  ind = find((t>tstart)&(t<tend));
  t = t(ind);
  amat = amat(ind,:);
end

plot(amat,t)

