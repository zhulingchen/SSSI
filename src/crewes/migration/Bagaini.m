function out=Bagaini(bin,Lin,vel)
% BAGAINI: find a piecewise constant velocity vel appropriate for PSPI
%
% out=Bagaini(bin,Lin,vel)
%
% See Bagaini's paper on choosing reference velocities for PSPI, either in
% Geophysics or conference abstracts.
%
% bin ... usually make this = size(vel,2)-1
% Lin ... not sure what this is, try 10
% vel ... input velocity model
%
% Code by R.J. Ferguson, not commented and opaque
% Cleaned up a bit by G.F. Margrave, still opaque
%
[rm,cm]=size(vel);
out=zeros(rm,cm);
for j=1:rm
	b=cm;
	L=Lin;
	for k=cm:-1:1
		if b>bin
			L=L-1;
			c=admit_vels(vel(j,:),L);
			P=prob_den(c,vel(j,:));
			S=stat_ent(P);
			b=binz(S,P);
		else
			k=1;
		end
	end
	v=bound_ref_vels(c,b,P);
	out(j,:)=piece_wise_itize(v,vel(j,:));
end

function c=admit_vels(v,L)
[rv,cv]=size(v(:));
cmin=min(v);
dv=(max(v)-cmin)/L;
c=zeros(1,L+1);
for k=1:L+1
	c(k)=cmin+(k-1)*dv;
end

function P=prob_den(c,v)
[rc,cc]=size(c(:));
L=rc-1;
P=zeros(1,L);
for k=1:L
	P(k)=length(find(and(c(k)<=v,c(k+1)<=v)));
end
P=P/sum(P);

function S=stat_ent(P)
S=-sum(P.*log(P));

function b=binz(S,P)
[rp,cp]=size(P(:));
b=max(1,min(rp,round(exp(S)+.5)));

function v=bound_ref_vels(c,b,P)
[rc,cc]=size(c(:));
y=zeros(1,rc);
y(1)=0;
y(rc)=1;
for j=2:rc-1
	y(j)=sum(P(1:j-1));
end
v=zeros(1,b+1);
v(1)=c(1);
for j=1:b
	for k=1:rc-1
		if and(y(k)<j/b,i/b<=y(k+1))
			v(j+1)=c(k)+(j/b-y(k))*(c(k+1)-c(k))/(y(k+1)-y(k));
		end
	end
end

function out=piece_wise_itize(v,vel)
[rv,cv]=size(v(:));
[rm,cm]=size(vel(:));
out=v(rv)+zeros(1,rm);
inds=0;
for j=1:rv-1
	inds=find(and(vel>=v(j),vel<=v(j+1)));
	out(inds)=mean(vel(inds));
end
