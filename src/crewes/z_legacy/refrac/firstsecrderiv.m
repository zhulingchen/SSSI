%first break vs offset
figure('menubar','none')
%shotradpickm=median(shotradpick);
%[b,a]=fir1(6,[0.1,0.2]);
%shotradpickf=filtfilt(b,a,shotradpick);
plot(shotoffset(1:189,:),shotradpick(1:189,:),'o');
xlabel('offset (m)')
ylabel('traveltime (ms)')
title('Refracted radial picks for shots #1-189')
%first derivative
[m,n]=size(shotoffset);
delx=shotoffset(:,1:n-40)-shotoffset(:,41:n);
rdelt=shotradpick(:,1:n-40)-shotradpick(:,41:n);  
avgx=(shotoffset(:,1:n-40)+shotoffset(:,41:n))/2;
rderiv=delx./rdelt;
absavgx=abs(avgx);
absrderiv=abs(rderiv);
figure('menubar','none')
plot(absavgx(1:189,:),absrderiv(1:189,:),'o');
axis([0,4000,0,2]);                     
xlabel('offset (m)')
ylabel('velocity (1000 m/s)')
title('First derivative of the radial picks for shots #1-189')
%second derivative
[m2,n2]=size(rderiv);
delx2=delx(:,1:n2-10)-delx(:,11:n2);
rdelt2=rdelt(:,1:n2-10)-rdelt(:,11:n2);
avgx2=(avgx(:,1:n2-10)+avgx(:,11:n2))/2;
rderiv2=delx2./rdelt2;
absavgx2=abs(avgx2);
%absrderiv2=abs(rderiv2);
figure('menubar','none')
plot(absavgx2(1:189,:),rderiv2(1:189,:),'o');
axis([0,4000,-0.2,0.2]);  
xlabel('offset (m)')
ylabel('acceleration (1000 m/s2)')
title('Second derivative of the radial picks for shots #1-189')
