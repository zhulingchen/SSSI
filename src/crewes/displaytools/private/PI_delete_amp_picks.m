function PI_delete_amp_picks(hfig)
%
% deletes the AMP_PICKS for the figure hfig
% Gets called either by PI_Close or PI_zoompick
%

global AMP_PICKS
nevents=length(AMP_PICKS);
knew=0;
new_amp_picks=cell(1,nevents);
for k=1:nevents
    pickstruc=AMP_PICKS{k};
    if(pickstruc.figurehandle~=hfig)
        knew=knew+1;
        new_amp_picks{knew}=AMP_PICKS{k};
    else
        %yes this really is necessary
        if(ishandle(pickstruc.handle))
            delete(pickstruc.handle);
        end
        if(ishandle(pickstruc.trajhandle))
            delete(pickstruc.trajhandle);
        end
        if(ishandle(pickstruc.texthandle))
            delete(pickstruc.texthandle);
        end
    end
end
AMP_PICKS=new_amp_picks(1:knew);
