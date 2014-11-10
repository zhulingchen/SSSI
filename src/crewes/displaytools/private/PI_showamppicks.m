function PI_showamppicks

global AMP_PICKS
h=get(gcf,'userdata');
hfig=gcf;
hamp_picks_lbl=h(18);

show=get(hamp_picks_lbl,'checked');
if(strcmp(show,'on'))
    show='off';
else
    show='on';
end
set(hamp_picks_lbl,'checked',show)

nevents=length(AMP_PICKS);

for k=1:nevents
    pickstruc=AMP_PICKS{k};
    if(pickstruc.figurehandle==hfig)
        set(pickstruc.texthandle,'visible',show);
    end
end