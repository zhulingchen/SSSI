function PI_PlotImageColorMap;
h=get(gcf,'userdata');
hmsg=h(4);
hi=h(5);
hscale=h(6);
hclip=h(7);

%get the data
seis=get(hi,'cdata');

%determine old scaling
dat=get(hscale,'userdata');
oldscaleopt=dat(1);
mxs=dat(2); mns=dat(3); smean=dat(4); stddev=dat(5);
colorview(gca,hi,mns,mxs,0)
plotimage('limboxrescale');
