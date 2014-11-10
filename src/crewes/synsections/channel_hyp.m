% CHANNEL: model a channel beneath a few layers
%
% low velocity channel beneath a v(z) medium
% Just run the script
dx=5; %cdp interval
dt=.002; %time sample rate
xmax=2500;tmax=1.0; %maximum line length and maximum time
x=0:dx:xmax; % x coordinate vector
t=0:dt:tmax; % t coordinate vector
v=3000; % velocity
z=v*t/2; % z coordinate vector

%initialize seismic matrix
seis=zeros(length(t),length(x));

%first event
z1=100;
seis=event_pwlinh(seis,t,x,v,[0 xmax],[z1 z1],[.1 .1]);
disp('first event of five done')

%second event
z2=200;
seis=event_pwlinh(seis,t,x,v,[0 xmax],[z2 z2],[-.1 -.1]);
disp('second event of five done')

%third event
z3=271;
seis=event_pwlinh(seis,t,x,v,[0 xmax],[z3 z3],[.05 .05]);
disp('third event of five done')

%fourth event
z4=398;
seis=event_pwlinh(seis,t,x,v,[0 xmax],[z4 z4],[-.04 -.04]);
disp('fourth event of five done')

%channel
width=100;
thk=50;
xchan=[xmax/2-width/2  xmax/2-.5*width/2 xmax/2+.5*width/2 xmax/2+width/2];
zchan=[z4 z4+thk z4+thk z4];
seis=event_pwlinh(seis,t,x,v,xchan,zchan,.04*ones(size(xchan)));
disp('fifth event of five done')

plotimage(seis,t,x);title('Unfiltered model')

%apply a filter

seis=sectfilt(seis,t,[10 5],[150 30]);

plotimage(seis,t,x);title('Filtered to 10-150 Hz.');
