function [plust] = avgplustime(fbcoord,td1,standard,plustreject,dev)
% Calculation of the average Plus Time at each receiver location
% Rejection of some of the Plus Time based on the standard deviation or
% on a constant difference from the average Plus Time values
recelev=refdata('get','recelev');
% receiver coordinates
coord=recelev(1,:);
% Define the average Plus Time matricies:
% coord: receiver coordinates
% plust: average Plus Time values
% std: standard deviation
% fold:  Plus Time fold
[tmp ncoords] = size(coord);
plust = NaN.*ones(3, ncoords);
foldplust = zeros(1, ncoords);
[nplust endhg] = size(td1);
allplust =NaN.*ones(nplust,ncoords);
avgplust=NaN.*ones(1,ncoords);
stdplust=NaN.*ones(1,ncoords);
% Order the Plus Time in term of coordinate location
for n=1:nplust
   fprintf(1,'Averaging for shots %d %d\n',td1(n,1), td1(n,2));
   i = td1(n,1);
   d = td1(n,3:endhg);
   dgoodi = find(~isnan(d));
   % Find the coordinate of the valid Plus Time in the shot i fbcoord  
   for val=dgoodi
	x = fbcoord(i,val);
	ind=find(coord == x);
	allplust(n,ind)=d(val);
   end
end
% Averaging of all the Plus Time according to their coordinate  
for n=1:ncoords
	validplust=find(~isnan(allplust(:,n)));
	[a b]= size(validplust);
    if( b ~=0 )
	avgplust(n)=mean(allplust(validplust,n));
	stdplust(n)=std(allplust(validplust,n));
	foldplust(n)=a;
    end
end
% Rejection of some Plus Time according to the standard deviation 
% or a constant limit
if (plustreject==1)
  for n=1:ncoords
	validplust=find(~isnan(allplust(:,n)));
	[a b]= size(validplust);
    if( b ~=0 )
	x=allplust(validplust,n);
	d=abs(x-avgplust(n));
	if (standard==1)
	   f=dev * stdplust(n);
	else
	   f=dev;
	end
	
	badplust=find(d>f);
	[a b]=size(badplust);
	if (b ~=0)
	   x(badplust)=NaN*badplust;
	   allplust(validplust,n) = x;
	   validplust=find(~isnan(allplust(:,n)));
	   [a b]= size(validplust);
    	   if( b ~=0 )
		avgplust(n)=mean(allplust(validplust,n));
		stdplust(n)=std(allplust(validplust,n));
		foldplust(n)=a;
	   end
	end
    end
  end
end
% Interpolation of the average Plus Time value to each receiver location
validplust=find(~isnan(avgplust));
avgplust=interpextrap(coord(validplust),avgplust(validplust),coord);
% Store all the information concerning the Plus Time values in a matrix of 4 rows
% first one: receiver coordinates
% second one: average Plus Time values
% third one: fold of the Plus Time values
% fourth one: standard deviation of the Plus time values
plust(1,:)=coord;
plust(2,:)=avgplust;
plust(3,:)=foldplust;
plust(4,:)=stdplust;
