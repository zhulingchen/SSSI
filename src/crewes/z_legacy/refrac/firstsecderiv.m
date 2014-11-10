%first break vs offset
figure('menubar','none')
%shotpickm=median(shotpick);
%[b,a]=fir1(6,[0.1,0.2]);
%shotpickf=filtfilt(b,a,shotpick);
plot(shotoffset(1:189,:),shotpick(1:189,:),'o');
xlabel('offset (m)')
ylabel('traveltime (ms)')
title('Refracted P-wave picks for shots #1-189')
%first derivative
[m,n]=size(shotoffset);
delx=shotoffset(:,1:n-50)-shotoffset(:,51:n);
delt=shotpick(:,1:n-50)-shotpick(:,51:n);  
avgx=(shotoffset(:,1:n-50)+shotoffset(:,51:n))/2;
deriv=delx./delt;
absavgx=abs(avgx);
absderiv=abs(deriv);
figure('menubar','none')
plot(absavgx(1:189,:),absderiv(1:189,:),'o');
axis([0,3500,2,4]);                     
xlabel('offset (m)')
ylabel('velocity (1000 m/s)')
title('First derivative of the P-wave  picks for shots #1-189')
%second derivative
[m2,n2]=size(deriv);
delx2=delx(:,1:n2-10)-delx(:,11:n2);
delt2=delt(:,1:n2-10)-delt(:,11:n2);
avgx2=(avgx(:,1:n2-10)+avgx(:,1:n2-10))/2;
deriv2=delx2./delt2;
absavgx2=abs(avgx2);
%absderiv2=abs(deriv2);
figure('menubar','none')
plot(absavgx2(1:189,:),deriv2(1:189,:),'o');
axis([0,4000,-0.5,0.5]);  
xlabel('offset (m)')
ylabel('acceleration (1000 m/s2)')
title('Second derivative of the P-wave picks for shots #1-189')
