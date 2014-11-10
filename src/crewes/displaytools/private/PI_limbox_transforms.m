function PI_limbox_transforms(hObject, eventdata, handles)
h=get(gcf,'userdata');
hmsg=h(2);
hi=h(5);
hlimbox=h(14);
nm=get(gcbo,'label');
httl=get(gca,'title');  ttl=get(httl,'string');
seis=get(hi,'cdata');
t=get(hi,'ydata');
x=get(hi,'xdata');
limfig=get(hlimbox,'userdata');
limdat=get(limfig{3},'userdata');
limmat=[];
for ii=1:4
    limmat=[limmat str2num(get(limdat(ii),'string'))];
end
imxdat=get(hi,'xdata'); imydat=get(hi,'ydata');
xdat=find(imxdat>=limmat(3) & imxdat<=limmat(4));
ydat=find(imydat>=limmat(1) & imydat<=limmat(2));
if(isempty(xdat))
    return
end
dddd=[xdat(1) xdat(end) ydat(1) ydat(end)];
seis=seis(ydat,xdat);
xdat=imxdat(xdat);
ydat=imydat(ydat);
if(isempty(seis))
    return
end
switch nm
    case 'F-K Amp Spectrum'
        stringinfo1='Generating F-K Amplitude Spectrum';
        stringinfo2='F-K Amplitude Spectrum';
        % limiting data that will FK transformed
        % [spec,f,kx]=fktran(seis,ydat,xdat,2^nextpow2(ydat),2^nextpow2(xdat),20);
        [spec,f,kx]=fktran(seis,ydat,xdat);
        xlbl='Wave Number';
        ylbl='Frequencey';
        file=['FK Spectrum of: ' ttl];
        init_image;
        h=get(gcf,'userdata');
        hscale=h(6);
        hclip=h(7);
        hmaster=h(10);
        for ii=1:length(h)
            if(ishandle(h(ii))&ii~=5)
                set(h(ii),'enable','on');
            end
        end
        newim=image(kx,f,abs(spec));
        colormap('gray');
        imagetype='Amp Spectrum';    % multiple different image types
    case 'F-X Amp Spectrum'
        stringinfo1='Generating F-X Amplitude Spectrum';
        stringinfo2='F-X Amplitude Spectrum';
        % [spec,f]= fftrl(seis,ydat,10,2);
        [spec,f]= fftrl(seis,ydat);
        xlbl='Distance';
        ylbl='Frequencey';
        file=['FX Spectrum of: ' ttl];
        init_image;
        h=get(gcf,'userdata');
        hscale=h(6);
        hclip=h(7);
        hmaster=h(10);
        for ii=1:length(h)
            if(ishandle(h(ii))&ii~=5)
                set(h(ii),'enable','on');
            end
        end
        spec=abs(spec);
        ky=(1:1:size(spec,2));
        newim=image(ky,f,spec);
        colormap('gray');
        imagetype='Amp Spectrum';    % multiple different image types
    case 'More to Come'
    case 'More to Come'
end
set(hmsg,'string',stringinfo1,'backgroundcolor',[1 1 1]);
hmsg=h(2);
set(hmsg,'string',stringinfo2,'backgroundcolor',[1 1 1]);
set(gca,'xaxislocation','bottom');
gtfigs=findobj(0,'type','figure','tag','PLOTIMAGEFIGURE');
nm=1;
xxnm=1;
adon='';
for ii=1:length(gtfigs)
    haxs=get(gtfigs(ii),'currentaxes');
    ttl=get(haxs,'title');
    dat=get(ttl,'userdata');
    if(~isempty(dat))
        xfile=dat{1};
        xnm=dat{2};
        if(strcmp(xfile,file));
            nm=nm+1;
            if(xnm>=nm)
                nm=xnm+1;
            end
            adon=['(' num2str(nm) ')'];
        end 
    end
end
title([file adon],'tag','PLOTIMAGE-TITLE','fontweight','bold',...
    'userdata',{file nm imagetype},'interpreter','none');
xlabel(xlbl,'horizontalalignment','right');
ylabel(ylbl,'tag','PLOTIMAGE-YLABEL');
h=get(gcf,'userdata');
h(5)=newim;
h(10)=hmaster;
set(hmaster,'visible','off');
mxs=full(max(max(abs(spec))));
%determine clipping
smean=full(mean(mean(spec)));
stddev=full(sqrt( sum(sum((spec-smean).^2 ) )...
    /prod(size(spec))));
clip=4;
smean2=smean;
mxs2=mxs;
stddev2=stddev;
scaleopt=1;
if(~isnan(clip))
    mxsprime=min([smean2+clip*stddev2,mxs2]);
end
mns=-mxsprime;
mns2=mns;
cm=get(gcf,'colormap');
[nkols,m]=size(cm);
seis = (spec-mns)/(mxsprime-mns)*(nkols-1)+1;
clear smat
ampflag=2;
set(h(6),'userdata',[scaleopt mxs2 mns2 smean2 stddev2]);
set(h(10),'userdata',[mxs smean stddev]);
set(gcf,'userdata',h);
