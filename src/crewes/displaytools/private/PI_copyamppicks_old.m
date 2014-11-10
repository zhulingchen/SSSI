function PI_copyamppicks(action)

global AMP_PICKS
persistent stuff

if(nargin<1)
    action='init';
end

if(strcmp(action,'init'))
    h=get(gcf,'userdata');
    hzoompick=h(9);
    set(hzoompick,'value',4);
    PI_zoompick;
    %find other plotimage windows with picks
    nevents=length(AMP_PICKS);
    hthisfig=gcf;%the figure window that the new picks will be in
    hotherfigs=zeros(1,20);%handles of the other plotimage windows
    eventfigs=zeros(1,100);%figure number for each event
    notherfigs=0;%number of other plotimage windows
    notherevents=0;%number of events in all other windows
    eventnames='';
    for k=1:nevents
        pickstruc=AMP_PICKS{k};
        hfig=pickstruc.figurehandle;
        if(hfig~=hthisfig)
            notherevents=notherevents+1;
            ind=find(hotherfigs,hfig);
            if(isempty(ind))
                notherfigs=notherfigs+1;
                hotherfigs(notherfigs)=hfig;
            end
            if(notherevents>1)
                eventnames=char(eventnames,[pickstruc.eventname ' in figure ' int2str(hfig)]);
            else
                eventnames=char([pickstruc.eventname ' in figure ' int2str(hfig)]);
            end
            eventfigs(notherevents)=hfig;
        end
    end
    if(isempty(eventnames))
        msgbox('No picked events in other Plotimage windows found');
        return
    end
    %clean up
    ind=find(eventfigs==0);
    eventfigs(ind)=[];
    
    %sort the event names by figure number
    [~,ind]=sort(eventfigs);
    eventnames=eventnames(ind,:);
    eventnames_rowvec=char(32*ones(1,numel(eventnames)+notherevents));
    kbegin=1;
    for k=1:notherevents
        name=[deblank(eventnames(k,:)) '|'];
        kend=kbegin+length(name)-1;
        eventnames_rowvec(kbegin:kend)=name;
        kbegin=kend+1;
    end
    %put up an askthings dialog
    q=char('Choose which event to copy','Choose picking option');
    a=char(eventnames_rowvec,...
        ['Pick at exact times of other event|Pick new event according to ' ...
        'specs of other event|Respecify picking using other event''s trajectory']);
    askthingsinit('plotimage(''copyamppicks2'')',q,a,[1 2],'Chose event and picking scheme');
    set(gcf,'windowstyle','modal');
elseif(strcmp(action,'pick')||strcmp(action,'pick2'))
    if(strcmp(action,'pick'))
        a=askthingsfini;
        if(a==-1)
            return;
        end
        neweventname=deblank(a(1,:));
        eventname=neweventname;
        ind=strfind(eventname,' in figure ');
        hfig=str2num(eventname(ind(1)+11:end));%figure number of event
        eventname=eventname(1:ind(1)-1);%name of event
        method=deblank(a(2,:));%picking method
        if(strcmp(method,'Respecify picking using other event''s trajectory'))
            %build the picker dialog as in ipick
            q=str2mat('Half-width of picking fairway (sec)',...
                'Pick what?');
            a=str2mat('0.020',...
              ['pick max(abs) amp|pick nearest peak of Hilbert env|pick nearest peak'...
               '|pick nearest trough|pick nearest +to- zero crossing'...
               '|pick nearest -to+ zero crossing|pick nearest zero crossing']);
            flags=[1 2];
            transferfcn='plotimage(''copyamppicks3'')';
            stuff.name=neweventname;
            stuff.hfig=hfig;
            askthingsinit(transferfcn,q,a,flags,['Picking Dialog for ' neweventname]);
            set(gcf,'windowstyle','modal');

            return
        end
    elseif(strcmp(action,'pick2'))
        a=askthingsfini;
        if(a==-1)
            return;
        end
        neweventname=stuff.name;
        eventname=neweventname;
        ind=strfind(eventname,' in figure ');
        eventname=eventname(1:ind(1)-1);%name of event
        hfig=stuff.hfig;
        delt=str2double(a(1,:));
        pickmode=a(2,:);
        method='3';
    end
    %ok now get pickstruc
    nevents=length(AMP_PICKS);
    for k=1:nevents
        pickstruc=AMP_PICKS{k};
        hfignow=pickstruc.figurehandle;
        if(hfig==hfignow)
            if(strcmp(eventname,pickstruc.eventname))
                break;
            end
        end
    end
    %turn pickmode into a flag
    if(~strcmp(method,'3'))
        pickmode=pickstruc.picktype;
        if(strcmp(pickmode,'pick max(abs) amp'))
            flag=1;
        elseif(strcmp(pickmode,'pick nearest peak of Hilbert env'))
            flag=2;
        elseif(strcmp(pickmode,'pick nearest peak'))
            flag=3;
        elseif(strcmp(pickmode,'pick nearest trough'))
            flag=4;
        elseif(strcmp(pickmode,'pick nearest +to- zero crossing'))
            flag=5;
        elseif(strcmp(pickmode,'pick nearest -to+ zero crossing'))
            flag=6;
        elseif(strcmp(pickmode,'pick nearest zero crossing'))
            flag=7;
        else
            error('logic failure in ipick, you''re screwed')
        end
        %turn method into a flag
        if(strcmp(method,'Pick at exact times of other event'))
            methodflag=1;
        elseif(strcmp(method,'Pick new event according to specs of other event'));
            methodflag=2;
        else
            error('unable to resolve method, you''re screwed')
        end
        if(methodflag==1 && (flag==5 || flag==6 || flag==7))
            methodflag=2;
        end
    else
        if(strcmp(pickmode,'pick max(abs) amp'))
            flag=1;
        elseif(strcmp(pickmode,'pick nearest peak of Hilbert env'))
            flag=2;
        elseif(strcmp(pickmode,'pick nearest peak'))
            flag=3;
        elseif(strcmp(pickmode,'pick nearest trough'))
            flag=4;
        elseif(strcmp(pickmode,'pick nearest +to- zero crossing'))
            flag=5;
        elseif(strcmp(pickmode,'pick nearest -to+ zero crossing'))
            flag=6;
        elseif(strcmp(pickmode,'pick nearest zero crossing'))
            flag=7;
        else
            error('logic failure in ipick, you''re screwed')
        end
        methodflag=3;
    end
    [seis,t,x]=plotimage_getseismic;%get the seismic
    xe=get(pickstruc.trajhandle,'xdata');
    te=get(pickstruc.trajhandle,'ydata');
    kol=get(pickstruc.trajhandle,'color');
    htraj=line(xe,te,'marker','*','color',kol,'linestyle','none');%plot trajectory
    if(methodflag==2);
        delt=pickstruc.delt;%fairway width
    elseif(methodflag==1)
        delt=0;
    end
    [ap,ae,tp,xp]=picker(seis,t,x,te,xe,delt,flag);
    hcntx=uicontextmenu;
    uimenu(hcntx,'label',deblank(neweventname));
    % Make a new pickstructure
    pickspec.figurehandle=gcf;
    pickspec.handle=line(xp,tp,'marker','.','linestyle','none','color',kol,...
        'uicontextmenu',hcntx);
    pickspec.eventname=deblank(neweventname);
    pickspec.picktype=pickmode;
    pickspec.trajhandle=htraj;
    pickspec.amppick=ap;
    pickspec.ampevent=ae;
    pickspec.tpick=tp;
    pickspec.xpick=xp;
    pickspec.delt=delt;
    
    %return to plotimage
    set(gca,'userdata',pickspec);
    plotimage('fromipick');
end
            
            
                