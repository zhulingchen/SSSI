function PI_DeletePickLine();
h=get(gcf,'userdata');
hmsg=h(2);
global PICKS 
delete(findobj(gcf,'type','line','tag','PICKMARKER'));
delete(findobj(gcf,'type','text','tag','PICKTEXT'));
hobj=gco;
for ii=1:size(PICKS,1)
    CheckFigure=PICKS{ii,1};
    if(CheckFigure==gcf)
        PicksPositions=PICKS{ii,2};
        PicksHandles=PICKS{ii,3};
        for jj=1:length(PicksHandles)
            if(PicksHandles(jj)==hobj)
                PicksPositions(jj,:)=[];
                delete(hobj);
                PicksHandles(jj)=[];
                PICKS{ii,2}=PicksPositions;
                PICKS{ii,3}=PicksHandles;
                break
            end
        end
    end
end
stringinfo=['Pick Line has been deleted'];
set(hmsg,'string',stringinfo,'backgroundcolor',[1 1 1]);
