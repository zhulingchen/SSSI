function unhideui(figno)
% UNHIDEUI ... restore user interface controls
%
% Makes visible all user interface controls in a figure
%
% G.F. Margave, CREWES, 2000
%

if(nargin<1) handin=gcf; end

flag=get(handin(1),'type');

if(strcmp(flag,'figure'))
   hhi=get(handin(1),'children');
   hh=nan*hhi;
   for kk=1:length(hhi)
	   tp=get(hhi(kk),'type');
       if(strcmp(tp,'uicontrol'))
           flg=get(hhi(kk),'visible');
           if(strcmp(flg,'off'))
               hh(kk)=hhi(kk);
               set(hhi(kk),'visible','on');
           end
       end
       if(strcmp(tp,'buttongroup'))
           flg=get(hhi(kk),'visible');
           if(strcmp(flg,'off'))
               hh(kk)=hhi(kk);
               set(hhi(kk),'visible','on');
               hk=get(hhi(kk),'children');
               for k=1:length(hk)
                   set(hk(k),'visible','on');
               end
           end
       end
       if(strcmp(tp,'uipanel'))
           flg=get(hhi(kk),'visible');
           if(strcmp(flg,'off'))
               hh(kk)=hhi(kk);
               set(hhi(kk),'visible','on');
               hk=get(hhi(kk),'children');
               for k=1:length(hk)
                   set(hk(k),'visible','on');
               end
           end  
       end
   end

   ind=find(isnan(hh));
   hh(ind)=[];
else
   for k=1:length(handin)
      set(handin(k),'visible','on');
   end
   hh=[];
end
