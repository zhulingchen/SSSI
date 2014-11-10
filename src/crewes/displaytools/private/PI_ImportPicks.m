function PI_ImportPicks
% Import picks from other figures in familly
global PICKS PICKCOLOR 
h=get(gcf,'userdata');
hmsg=h(2);
% getting name of figure pulling
ttl=get(gcbo,'label');
xx=findstr(ttl,':');
ttl=ttl(xx(1)+2:end);
NewPicksData=get(gcbo,'userdata');
MasterFigure=gcf;
CheckPicks=findobj(gcf,'type','line','tag','PICKS');
if(~isempty(CheckPicks))
    CheckwUser=questdlg('Picks already present.',...
        'Picks Present Alert','Delete Present Picks','Merge Picks','Cancel','Cancel');
    switch CheckwUser
    case 'Delete Present Picks'
        for ii=1:size(PICKS,1)
            if(PICKS{ii,1}==MasterFigure)
                delete(findobj(gcf,'type','line','tag','PICKS'));
                PICKS{ii,2}=[];
                PICKS{ii,3}=[];
                break
            end
        end
    case 'Merge Picks'
    case 'Cancel'
        return
    end
end
xcheck=get(gca,'xlim');
ycheck=get(gca,'ylim');
for ii=1:size(PICKS,1)
    if(PICKS{ii,1}==MasterFigure)
        PicksData=PICKS{ii,2};
        PicksHandles=PICKS{ii,3};
        for jj=1:size(NewPicksData,1)
            xdatsrt=sort([NewPicksData(jj,1) NewPicksData(jj,3)]);
            ydatsrt=sort([NewPicksData(jj,2) NewPicksData(jj,4)]);
            if(xcheck(1)>=xdatsrt(1)|xcheck(2)<=xdatsrt(2)|ycheck(1)>=ydatsrt(1)|ycheck(2)<=ydatsrt(2))
            else
                hpick=line([NewPicksData(jj,1) NewPicksData(jj,3)],...
                    [NewPicksData(jj,2) NewPicksData(jj,4)],[1 1],'linewidth',2,...
                    'color',PICKCOLOR,'buttondownfcn','plotimage(''picklinemenu'')',...
                    'userdata','','tag','PICKS');
                PicksData=[PicksData; NewPicksData(jj,1) NewPicksData(jj,2) NewPicksData(jj,3) NewPicksData(jj,4)];
                PicksHandles=[PicksHandles;hpick];
            end
        end
        PICKS{ii,2}=PicksData;
        PICKS{ii,3}=PicksHandles;
        break
    end
end
xx=length(PICKS{ii,3});
delete(findobj(gcf,'type','line','tag','PICKMARKER'));
delete(findobj(gcf,'type','text','tag','PICKTEXT'));
stringinfo=[num2str(xx) ' Picks have been imported from "' ttl '"'];
set(hmsg,'string',stringinfo,'backgroundcolor',[1 1 1]);
