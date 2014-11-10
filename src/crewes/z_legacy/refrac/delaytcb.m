function delaytcb(action)
if( nargin < 1 )
   action = 'init';
end
if( strcmp(action,'init'))
   q=str2mat('Delay time analysis according to ?',...
	     'Enter the coordinate limit ?',...  
	     'Use a constant velocity for the second layer?',...
             'Enter the constant velocity (m/s):',...
	     'Autorejection of Delay times ?',...
	     'Autorejection based on constant or standard devation limit?',...
	     'Enter limit:');
   a=str2mat('CVPavg|Cte','250 7750','No|Yes','1600','Yes|No','Standard deviation|Constant limit','1')
   askthingsinit('delaytcb(''answer'')',q,a,[1 0 1 0 1 0 0],...
		 'Parameter for the Delay time analysis');
elseif( strcmp(action,'answer'))
   a=askthingsfini;
   [strings tmp] = size(a);
   if(strcmp( deblank(a(1,:)),'CVPavg'))
	cvpavg=refdata('get','cvpavg');
	ind1=find(~isnan(cvpavg(:,2)));
	s=ind1(1);
	slim=cvpavg(s,2);
	
	ind2=find(~isnan(cvpavg(:,1)));
	e=ind2(length(ind2));
	elim=cvpavg(e,1);
    else
	lim = sscanf(a(2,:), '%d %d');
	slim=lim(1);
	elim=lim(2);
    end
   if(strcmp( deblank(a(3,:)),'No'))
	v2rec=refdata('get','v2rec');
   else
	recelev=refdata('get','recelev');
	v2 = str2num(a(4,:));
	v2 = v2/1000;
        n = length(recelev);
        v2rec = v2*ones(1,n);
   end
   if(strcmp( deblank(a(5,:)),'Yes'))
	delayreject=1;
   else
	delayreject=0;
   end
   if(strcmp( deblank(a(6,:)),'standard deviation'))
	standard=1;
   else
	standard=0;
   end
   dev = str2num(a(7,:));
   refdata('set','standard',standard);
% Function calling the Delay time analysis function
  fbtime=refdata('get','fbtime');
  fbcoord=refdata('get','fbcoord');
  shotcoord=refdata('get','shotcoord');
  cvpavg=refdata('get','cvpavg');
  nshots=refdata('get','nshots');
  recelev=refdata('get','recelev');
  plust=refdata('get','plust');
delay = delayt(fbtime,fbcoord,cvpavg,v2rec,shotcoord,nshots,recelev,slim,elim,plust);
  d=length(recelev);
  for n=1:d
    good=find(~isnan(delay(:,n)));
    if (length(good)>0)
     plust(2,n)=2*mean(delay(good,n));
     plust(3,n)=length(good);
     plust(4,n)=2*std(delay(good,n));
     if (delayreject==1)
       if (standard==1)
	 f=dev * (plust(4,n))/2;
       else
	 f=dev;
       end
       d=abs(delay(good,n)-(plust(2,n))/2);
       badplust=find(d>f);
       [a b]=size(badplust);
	if (b ~=0)
	   delay(good(badplust),n)=NaN*badplust;
	   good=find(~isnan(delay(:,n))); 	  
    	   if (length(good)>0)
             plust(2,n)=2*mean(delay(good,n));
             plust(3,n)=length(good);
             plust(4,n)=2*std(delay(good,n));
	   end
	end
     end
    end
  end
  refdata('set','plust',plust);
end
