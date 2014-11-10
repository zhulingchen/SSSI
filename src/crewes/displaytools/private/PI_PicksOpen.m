function PI_PicksOpen
global PICKS PICKCOLOR 
mainax=findobj(gcf,'type','axes','tag','MAINAXES');
posax=findobj(gcf,'type','axes','tag','POSITIONAXES');
set(gcf,'currentaxes',mainax);
MasterFigure=gcf;
CheckPicks=findobj(gcf,'type','line','tag','PICKS');
if(~isempty(CheckPicks))
    CheckwUser=questdlg('Picks already present.',...
        'Picks Present Alert','Delete Picks','Merge Picks','Cancel','Cancel');
    switch CheckwUser
    case 'Delete Picks'
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
[file,path]=myuifile(gcf,'.mat','Please Choose Picks File','get');
checkpicks=load([path file]);
nm=fieldnames(checkpicks);
if(length(nm)>=2)
    % if is the picks were previously saved in plot image, there should 
    % only be one file
    checkpicks=matinterogatorlite([path file],1);
    if(isempty(checkpicks))
        return
    end
end
if(isstruct(checkpicks))
    nm=fieldnames(checkpicks);
    NewPicksData=getfield(checkpicks,nm{1,:});
elseif(iscell(checkpicks))
    NewPicksData=checkpicks{:,:};
end
if(~(4==size(NewPicksData,2)))
    return
end
xcheck=get(mainax,'xlim');
ycheck=get(mainax,'ylim');
for ii=1:size(PICKS,1)
    if(PICKS{ii,1}==MasterFigure)
        PicksData=PICKS{ii,2};
        PicksHandles=PICKS{ii,3};
        for jj=1:size(NewPicksData,1)
            for kk=1:4
                sp=[];
                if(~isnumeric(NewPicksData(jj,kk)))
                    sp=1; % if any part of New Picks are not numberica, skip
                    break
                end
            end
            xdatsrt=sort([NewPicksData(jj,1) NewPicksData(jj,3)]);
            ydatsrt=sort([NewPicksData(jj,2) NewPicksData(jj,4)]);
            if(xcheck(1)>=xdatsrt(1)|xcheck(2)<=xdatsrt(2)|ycheck(1)>=ydatsrt(1)|ycheck(2)<=ydatsrt(2))
            else
                if(isempty(sp))
                    hpick=line([NewPicksData(jj,1) NewPicksData(jj,3)],...
                        [NewPicksData(jj,2) NewPicksData(jj,4)],[1 1],'linewidth',2,...
                        'color',PICKCOLOR,'buttondownfcn','plotimage(''picklinemenu'')',...
                        'userdata','','tag','PICKS');
                    PicksData=[PicksData; NewPicksData(jj,1) NewPicksData(jj,2) NewPicksData(jj,3) NewPicksData(jj,4)];
                    PicksHandles=[PicksHandles;hpick];
                end
            end
        end
        PICKS{ii,2}=PicksData;
        PICKS{ii,3}=PicksHandles;
        break
    end
end
delete(findobj(gcf,'type','line','tag','PICKMARKER'));
delete(findobj(gcf,'type','text','tag','PICKTEXT'));
return
