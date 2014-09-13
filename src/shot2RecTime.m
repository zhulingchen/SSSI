function it = shot2RecTime(travelTime,ixs,ixr,dt,nx)
if nx < size(travelTime,2)
    it = round(travelTime(:,1:ixs,ixs)/dt)+1 + round(travelTime(:,1:ixs,ixr)/dt)+1;
else
    it = round(travelTime(:,:,ixs)/dt)+1 + round(travelTime(:,:,ixr)/dt)+1;
end