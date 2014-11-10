function PI_ClosePlotImage
% Closing Plot image
global PICKS PICKCOLOR ZOOM_VALUE ZOOM_LOCKS
limfig=findobj(0,'type','figure','tag','MVLINESMEASUREMENTS');
h=get(gcf,'userdata');
hlimbox=h(14);
limdat=get(hlimbox,'userdata');
masterfigure=findobj(gcf,'type','uimenu','tag','PLOTIMAGEMASTERFIGURE');
mf=get(masterfigure,'parent');
mf=get(mf,'parent');
if(~isempty(mf))
    button = questdlg('Are you sure you want to quit?',...
        'Closing Plotimage','Yes','No','Yes');
    switch button
    case 'Yes'
    case 'No' 
        return
    end
end
if(size(PICKS,1)>=2)
    NewPICKS=cell(1,3);
    jj=1;
    for ii=1:size(PICKS,1)
        CheckFigure=PICKS{ii,1};
        if(~isempty(PICKS{ii,1}))
            if(CheckFigure==gcf)
            else
                NewPICKS{jj,1}=PICKS{ii,1};
                NewPICKS{jj,2}=PICKS{ii,2};
                NewPICKS{jj,3}=PICKS{ii,3};
                jj=jj+1;
            end
        end
    end
    PICKS=NewPICKS;
else
    PICKS=[];
end
if(~isempty(ZOOM_LOCKS))
    NewZOOM_LOCKS=[];
    for ii=1:size(ZOOM_LOCKS,1)
        if(~isempty(find(ZOOM_LOCKS(ii,:)==gcf)))
        else
            NewZOOM_LOCKS=[NewZOOM_LOCKS,ZOOM_LOCKS(ii,:)];
        end
    end
    ZOOM_LOCKS=NewZOOM_LOCKS;
end
if(~isempty(limdat))
    if(ishandle(limdat{3}))
        delete(limdat{3});
    end
end
% if(~isempty(masterfigure))
%     % if it is masterfigure, deleteing the other plotimages spawned
%     delete(findobj(0,'type','figure','tag','PLOTIMAGEFIGURE'));
%     clear global;
%     return
% end
delete(gcf);
% if(isempty(findobj(0,'type','figure','tag','PLOTIMAGEFIGURE')));
%     clear global;
% end
