function editdepth(edit,win)
% Edit the layer1 depth by draging layer1-layer2 interface 
% (values at each receivers) 
% Optional median filter applicable to depth values
depth=refdata('get','depth');
recelev=refdata('get','recelev');
set(gcf,'units','pixels','position',[0 0 864 576]);
% Median filter
if (edit==0)
  depth(2,:)=medfilt1(depth(2,:),win);
end
depthelev=recelev(2,:)-depth(2,:);
hold on;
d1handle = plot(depth(1,:),recelev(2,:),'color','y','linestyle','-')
d2handle = plot(depth(1,:),depthelev,'color','g','linestyle','-')
ylabel('Elevation (m)')
title('Depth model (surface=yellow; layer 1-2 interface=green)')
xlabel('Coordinate (m)')
% Call the edit function
editveldepth('init',[d2handle], 'depth');
