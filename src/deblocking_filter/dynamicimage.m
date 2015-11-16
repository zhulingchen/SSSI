function ttbin=dynamicimage(stackeddataM2int,nrpixels,imag_col)

if nargin==1
    %nrpixels=16;
    nrpixels=30;
    %nrpixels = 64;
    %nrpixels=128;
    %imag_col=255;
    imag_col = 255;
end
if nargin==2
    imag_col=64;
end

sc=size(stackeddataM2int,2);
begin_im=zeros(1,sc);end_im=zeros(1,sc);
for s=1:sc
    begin_im(s)=find(isnan(stackeddataM2int(:,s))==0,1,'first');
    end_im(s)=find(isnan(stackeddataM2int(:,s))==0,1,'last');
end
if isempty(begin_im);startdepthbin=1;end
if isempty(end_im);enddepthbin=size(stackeddataM2int,1);end

startdepthbin=min(begin_im);
enddepthbin=max(end_im);

W=nrpixels;
Ltot=enddepthbin-startdepthbin+1;
j=0;clear lvlsM2 pixnr;pixnr=zeros(1,floor(Ltot/W-1));
for i=startdepthbin:W:enddepthbin
    j=j+1;E=min(i+W-1,enddepthbin);pixnr(j)=round(i+W/2);
    lvlsM2dyn(1:imag_col,j)=colorbinlevels(stackeddataM2int(i:E,1:sc),imag_col,sc);
end
y=zeros(imag_col,length(stackeddataM2int));
for i=1:imag_col
    y(i,1:pixnr(1))=lvlsM2dyn(i,1);
    y(i,pixnr(end):enddepthbin)=lvlsM2dyn(i,end);
    if length(pixnr)>1
        y(i,pixnr(1):pixnr(end))=interp1(pixnr,lvlsM2dyn(i,:),pixnr(1):pixnr(end),'linear');
    end
end

x=stackeddataM2int;
clrs=y;

% %%%%%%%%%%%%%%make 64 step histogram data for coloring
startdpt=1;
[enddpt dummy]=size(x);
[nrclrs dummy]=size(clrs);
ttbin=zeros(enddpt-startdpt+1,sc,'int16');

for i=1:sc
    for j=startdpt:enddpt
      for k=nrclrs-1:-1:1
          if x(j,i)>clrs(k,j) %-1e-3 %if I use x=>clrs, it does not work properly 
              ttbin(j,i)=k;
              break
          end
      end
    end
end

[a b]=find(x==Inf);
for i=1:length(a);ttbin(a(i),b(i))=nrclrs;end
%[a b]=find(x==0);
%for i=1:length(a);ttbin(a(i),b(i))=nrclrs;end
[a2 b2]=find(isnan(x));
for i=1:length(a2);ttbin(a2(i),b2(i))=nrclrs;end
%if length(a)==0 & length(a2)==0;ttbin(end,end)=nrclrs;end

end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ttbinlevels=colorbinlevels(x,ncolors,sc)

% %%%%%%%%%%%%%%make step histogram data for coloring
startdpt=1;
[enddpt dummy]=size(x);
dptlong=zeros(sc*(enddpt-startdpt+1)-length(find(x==Inf))-length(find(isnan(x)))-length(find(x==0)),1);
m=0;
% drop all zeros and Nan and Inf's from the array
for i=1:sc
    for j=startdpt:enddpt
        if x(j,i)~=Inf & isnan(x(j,i))~=1 & x(j,i)~=0;
        m=m+1;
        dptlong(m)=x(j,i);
        end
    end
end
sizel=m;
i=m;
remove_val=zeros(1,m);l=0;
while i>1
    if dptlong(i-1)==dptlong(i)
        l=l+1;
        remove_val(l)=i;
    end
    i=i-1;
end
remove_val(l+1:end)=[];
dptlong(remove_val)=[];
sizel=length(dptlong);

%sort the resulting values and find the values where the color needs to
%change
% if sizel~=0
if sizel>ncolors
    sortdptlong=sortrows(dptlong,1);
    for i=1:ncolors
    clrs(i)=sortdptlong(round(i*sizel/ncolors));
    end
    ttbinlevels=clrs;
else
    ttbinlevels(1:ncolors)=-1;
end

end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
