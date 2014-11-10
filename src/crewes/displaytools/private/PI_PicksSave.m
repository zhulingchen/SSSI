function PI_PicksSave
global PICKS
hfig=gcf;
for ii=1:size(PICKS,1)
    if(PICKS{ii,1}==hfig)
        picksdata=PICKS{ii,2};
        break
    end
end
if(isempty(picksdata))
    return
end
[file,path]=myuifile(gcf,'.mat','Save Pick File','put');
if(isempty(file))
    return
end
save([path file],'picksdata');
