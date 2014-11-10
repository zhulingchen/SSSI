% Marmousi: finite difference modelling of Marmousi
%
% Marmousi finite difference modelling
% Just run the script

load marmousi_dz10
%load marmousi_mod
dx=x(2)-x(1);; %cdp interval
xmax=max(x);zmax=max(z); %maximum line length and maximum depth
vlow=min(min(vel));
vhigh=max(max(vel));

%plot the velocity model
plotimage(vel-mean(mean(vel)),z,x)
xlabel('meters');ylabel('meters')
velfig=gcf;

disp(' Type 1 for a shot record, 2 for a VSP, 3 for exploding reflector, <CR> to end')
r=input(' response -> ')

while ~isempty(r)
switch(r)
case{1}
	%do a shot record
	disp(' Type 1 for 2nd order Laplacian, 2 for 4th order Laplacian')
	lap=input(' response-> ');
	dt=.004; %temporal sample rate
	dtstep=.0005;%time step size
	tmax=2*zmax/vlow; %maximum time
	snap1=zeros(size(vel));
	xshot=max(x)/3;
	snap2=snap1; ix=near(x,xshot);
	snap2(1,ix(1))=1;
	[shotf,shot,t]=afd_shotrec(dx,dtstep,dt,tmax, ...
 		vel,snap1,snap2,x,zeros(size(x)),[5 10 40 50],0,lap);
	figure(velfig)
	line(xshot,0,1,'marker','*','markersize',6,'color','r');
	text(xshot+.02*(max(x)-min(x)),0,1,'shot point','color','r');
	%plot the seismogram
	plotimage(shotf,t,x-xshot)
   xlabel('offset');ylabel('seconds');
case{2}
	%do a vsp reflector model
	disp(' Type 1 for 2nd order Laplacian, 2 for 4th order Laplacian')
	lap=input(' response-> ');
	dt=.004; %temporal sample rate
	dtstep=.0005;%time step size
	tmax=2*zmax/vlow; %maximum time
	snap1=zeros(size(vel));
	xshot=max(x)/4;
	snap2=snap1; ix=near(x,xshot);
	snap2(1,ix(1))=1;
	izrec=near(z,max(z)/4,3*max(z)/4);
	zrec=z(izrec);
	xrec=(max(x)/3)*ones(size(zrec));
	[vspf,vsp,t]=afd_shotrec(dx,dtstep,dt,tmax, ...
 		vel,snap1,snap2,xrec,zrec,[5 10 40 50],0,lap);
	figure(velfig)
	line(xshot,0,1,'marker','*','markersize',6,'color','g');
	text(xshot+.02*(max(x)-min(x)),0,1,'vsp shot','color','g');
	line(xrec,zrec,ones(size(zrec)),'marker','v','markersize',6,'color','g');
	%plot the vsp
	plotimage(vspf',zrec,t)%vsp's are normally plotted with depth as vert coordinate
   ylabel('depth');xlabel('seconds');
case{3}
	%do an exploding reflector model
	disp(' Type 1 for 2nd order Laplacian, 2 for 4th order Laplacian')
	lap=input(' response-> ');
	dt=.004; %temporal sample rate
	dtstep=.001;
	tmax=2*zmax/vlow; %maximum time
	[seisfilt,seis,t]=afd_explode(dx,dtstep,dt,tmax, ...
 		vel,x,zeros(size(x)),[5 10 40 50],0,lap);
	%plot the seismogram
	plotimage(seisfilt,t,x)
   xlabel('meters');ylabel('seconds')
otherwise
	disp('invalid selection')
end

disp(' Type 1 for a shot record, 2 for a VSP, 3 for exploding reflector, <CR> to end')
r=input(' response ->')

end


