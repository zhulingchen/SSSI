function hh=hideui(handin)
% HIDEUI ... hide (or restore) user interface controls
%
% hh=hideui(figno) ... hides the uicontrols on figure figno
% hideui(hh) ... restores the uicontrols whos handles are in hh
%
% G.F. Margave, CREWES, 2000
%

if(nargin<1); handin=gcf; end

flag=get(handin(1),'type');

if(strcmp(flag,'figure'))
   hhi=get(handin(1),'children');
   hh=nan*ones(size(hhi));
   for kk=1:length(hhi)
	   tp=get(hhi(kk),'type');
	   if(strcmp(tp,'uicontrol'))
		   flg=get(hhi(kk),'visible');
            if(strcmp(flg,'on'))
             hh(kk)=hhi(kk);
             set(hhi(kk),'visible','off');
            end
       end
       if(strcmp(tp,'buttongroup'))
           flg=get(hhi(kk),'visible');
           if(strcmp(flg,'on'))
               hh(kk)=hhi(kk);
               set(hhi(kk),'visible','off');
               hk=get(hhi(kk),'children');
               for k=1:length(hk)
                   set(hk(k),'visible','off');
               end
           end
       end
       if(strcmp(tp,'uipanel'))
           flg=get(hhi(kk),'visible');
           if(strcmp(flg,'on'))
               hh(kk)=hhi(kk);
               set(hhi(kk),'visible','off');
               hk=get(hhi(kk),'children');
               for k=1:length(hk)
                   set(hk(k),'visible','off');
               end
           end

       end
           
       if(strcmp(tp,'axes'))
          axflag=get(hhi(kk),'tag');
          if(strcmp(axflag,'POSITIONAXES'))  
               hhc=get(hhi(kk),'children');
               set(hhi(kk),'visible','off');
               hh(kk)=hhi(kk);
               for kkk=1:length(hhc)
                   if(strcmp(get(hhc(kkk),'type'),'image'))
                       flg=get(hhc(kkk),'visible');
                       if(strcmp(flg,'on'))
                         set(hhc(kkk),'visible','off');
                       end
                   end
               end
          end
          if(strcmp(axflag,'MAINAXES'))
              psn=get(hhi(kk),'position');
              if(psn(1)==.2910)
                  psn(1)=0.100;
                  psn(3)=.8440;
              end
              set(hhi(kk),'position',psn);
              hh(kk)=hhi(kk);
          end
       end
   end

   ind=find(isnan(hh));
   hh(ind)=[];
else
   for k=1:length(handin)
      tg=get(handin(k),'tag');
      if(strcmp(tg,'MAINAXES'))
          psn=get(handin(k),'position');
              psn(1)=0.291;
              psn(3)=.6860;
          set(handin(k),'position',psn);
      elseif(strcmp(tg,'POSITIONAXES'))
          hhc=get(handin(k),'children');
           set(handin(k),'visible','on');
           for kkk=1:length(hhc)
               if(strcmp(get(hhc(kkk),'type'),'image'))
                   flg=get(hhc(kkk),'visible');
                   if(strcmp(flg,'off'))
                     set(hhc(kkk),'visible','on');
                   end
               end
           end
      else
          set(handin(k),'visible','on');
      end
   end
   hh=[];
end
	
