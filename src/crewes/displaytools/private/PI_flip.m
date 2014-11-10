function PI_flip();
h=get(gcf,'userdata');
hflip=h(1);
hmsg=h(4);
pol=get(hflip,'userdata');
pol=-1*pol;
clr=get(gcf,'colormap');
colormap(flipud(clr));
set(hflip,'userdata',pol);
if(pol==1)
    set(hmsg,'string','Polarity Normal');
elseif(pol==-1)
    set(hmsg,'string','Polarity Reversed');
end

