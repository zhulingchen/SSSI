% CHANGE_GRID_SPACING ... example script to interpolate a velocity model
% interpolate different grid spacing in velocity matrix

% grid spacing 20 m = change to 10m (halve)

load vmatrix;

x=1:1:300;
y=1:1:150;
y=y';

x1=0.5:0.5:300;
y1=0.5:0.5:150;
y1=y1';

vnewgrid = interp2(x,y,vmatrix,x1,y1);

% grid spacing 20m = change to 5 m (quarter)

load vmatrix;

x=1:1:300;
y=1:1:150;
y=y';

x2=0.25:0.25:300;
y2=0.25:0.25:150;
y2=y2';

vnewgrid2=interp2(x,y,vmatrix,x2,y2);

for j=1:3
    vnewgrid2(j,:)=vnewgrid2(4,:);
    vnewgrid2(:,j)=vnewgrid2(:,4);
end


