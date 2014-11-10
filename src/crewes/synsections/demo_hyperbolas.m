msg='Enter hyperbola spacing parameter as an integer between 1 and 50. <cr> to end: ';
ndelx=round(input(msg));
while(~isempty(ndelx))
	v=2000;dx=5;dt=.004;%basic model parameters
	x=0:dx:3000;%x axis
	t=0:dt:1.5;%t axis
	xcntr=max(x)/2;
	seis5=zeros(length(t),length(x));%allocate seismic matrix
	seis5=event_diph2(seis5,t,x,v,0,500,1000,ndelx,0,.1);
    disp('Working ...')
	seis5=event_diph2(seis5,t,x,v,500,xcntr-500,1000,ndelx,-45,.1);
    disp('Working ...')
	seis5=event_diph2(seis5,t,x,v,xcntr-500,xcntr+500,500,ndelx,0,.1);
    disp('Working ...')
	seis5=event_diph2(seis5,t,x,v,xcntr+500,max(x)-500,500,ndelx,45,.1);
    disp('Still at it ...')
	seis5=event_diph2(seis5,t,x,v,max(x)-500,max(x),1000,ndelx,0,.1);
    disp('Wrapping up ...')
	[w,tw]=ricker(dt,40,.2);%make ricker wavelet
	seis5=sectconv(seis5,t,w,tw);%apply wavelet
    plotimage(seis5,t,x);
    title(['Hyperbola spacing parameter ' int2str(ndelx)])
	ndelx=round(input(msg));
end
