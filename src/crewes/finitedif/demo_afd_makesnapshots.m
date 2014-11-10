% we will make a movie of a shot record in the marmousi model

%specify the major parameters
dx=10;%spatial grid size (m)
dt=.001;%temporal step size (s)
%Note: dx and dt must be chosen to satisfy the stability condition for the
%4th order (in space) second order (in time) scalar wave equation. This
%means that max(vel(:))*dt/dx<=sqrt(3/8) or the code will not run. 
tmax=2;%record length (s)
fdom=30;%desired dominant frequency (Hz)
dtmovie=.1;%the time sample rate of the movie (s)
laplacian=2;
boundary=2;

% make the velocity model
[vel,x,z]=marmousimodel(dx);%dx must be either 5, 10 or 25

%specify smapshots
xshot=mean(x);%x coordinate of the shot
zshot=0;%z coordinate of the shot
ixshot=near(x,xshot);
izshot=near(z,zshot);
snap1=zeros(size(vel));
snap2=snap1;
snap2(izshot,ixshot)=1;

%make a causal wavelet
[w,tw]=wavemin(dt,fdom,tmax);

%specify the snashop times
tsnaps=0:dtmovie:tmax;

%make the movie
snaps=afd_makesnapshots(dx,dt,vel,snap1,snap2,tsnaps,laplacian,boundary,w);

%make the snapshot titles
titles=cell(size(tsnaps));
for k=1:length(tsnaps)
    titles{k}=[' Snapshot at time ' num2str(tsnaps(k)) ' s'];
end

%show the results with plotgathers
plotgathers(snaps,x,z,'distance (m)','depth (m)',titles);

%show the results with plotsnaps
plotsnaps(snaps,vel,x,z,'distance (m)','depth (m)',titles,'Marmousi model');

%save the result for later
save marmsnaps snaps x z titles vel dx dt fdom w tw

%% Replot a saved result with a different colormap
load marmsnaps
plotsnaps(snaps,vel,x,z,'distance (m)','depth (m)',titles,'Marmousi model','jet');