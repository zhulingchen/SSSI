function [fbtime, fbcoord, shotcoord, nshots, nrecs, recelev, shotelev, uphole] = loadsynt
% Load files from synthetic data
% Data is loaded as five files:
%  fbpick:    (trace#) refracted arrival times
%  shotcoord: (shot#)  shot coordinate
%  reccoord:  (rec#)   receiver coordinate
%  shotelev:  (shot#)  shot elevation
%  recelev:   (rec#)   receiver elevation
 
more off
fprintf(1,'Loading offset\n');
% Load the coordinate file
load '/vol/matlab/toolbox/crewes/refrac/stnxcoordelev.txt';
% Assign array values from the file 
stn_vec = stnxcoordelev(:,1);
xcoord_vec = stnxcoordelev(:,2);
elev_vec = stnxcoordelev(:,3);
n=1:101;% shot numbers
i=96 + n*4;%shot station numbers
[a b]=size(xcoord_vec);
recelev=NaN * ones(2,a);
recelev(1,:)=xcoord_vec';
recelev(2,:)=elev_vec';
nshots=101;
nrecs=96;
fbtime = NaN * ones(nshots,nrecs);
station = zeros(nshots,nrecs);
fbcoord = zeros(nshots,nrecs);
uphole = zeros(nshots,1);
ngap=2;
% Open file for reading
fh = fopen('/vol/matlab/toolbox/crewes/refrac/fbpickshotrecs.txt','r');
if( fh == -1 ) 
   error('Could not open file');
end
% Fill the fbtime  matrix with the refracted arrival times 
% Fill the fbcoord matrix with corresponding receiver coordinates
% (rows=shot number: colums=consecutive trace number) 
% Fill the shotelev and shotcoord arrays
for n=1:nshots
   row = fscanf(fh, '%d %d %d %d %d %d %d %d',8);
   fprintf(1,'Reading %d\n',n);
   sin = row(1);
   nlive = row(2);
   startlive = row(3);
   i=96 + n*4;
   shotcoord(n)=xcoord_vec(4*n-3);
   shotelev(n)=elev_vec(4*n-3)';
   station(n,:) = [i-nrecs/2:i-ngap/2 i+ngap/2:i+nrecs/2];	
   goodones = find(station(n,:) > 99 & station(n,:) < 501);
   badones = find(station(n,:) < 100 | station(n,:) > 500);
   fbcoord(n,badones) = NaN * ones(size(badones));
   fbcoord(n,goodones) = xcoord_vec(station(n,goodones)-99)';
   for row=0:7
      data = fscanf(fh, '%d %d %d %d %d %d %d %d %d %d %d %d',12);
      fbtime(n,(1+row*12:(row+1)*12)) = data';
   end
end
shotcoord=shotcoord';
% Search for zeros in the fbtime matrix and replace them by NaNs.
for n=1:nshots
	ind=find(fbtime(n,:)==0);
	nind=size(ind);
	fbtime(n,ind)=NaN*ones(nind);
end
return;
