function [cvpi,cvpj,tdiff1, tdiff2, mddiff1, mddiff2, ...
          deltf1, deltf2, delts1, delts2, ...
          start1indexj, end1indexj, start2indexj, end2indexj] = ...
   autopick(gs1,gs2,shotcoord,fbcoord,fbtime,...
            window,nshots,nd,offsetrange1, offsetrange2, windmn)
% This function picks the cross over point coordinates according to 
% the maximum of the second derivative of the traveltime difference
% The cvpi correspond to the left cross over point of shot i
% The cvpj correspond to the right cross over point of shot j
cvpj=NaN .* ones(nshots,nshots);
cvpi=NaN .* ones(nshots,nshots);
% Sometimes the gs1 vector is a row, sometimes a column.  Handle this.
[a b] = size(gs1);
nshotpairs = max(a, b);
% Loop throught the specified shot pair(s)
for k=1:nshotpairs
   i = gs1(k);
   j = gs2(k);
   fprintf(1,'autopick: computing for shotpairs %d %d\n',i,j);
   % First call the traveltime substraction function to calculate
   % the traveltime difference on each side of the shot pair locations 
   
   tdiff1=[]; tdiff2=[]; start1indexj=[]; start2indexj=[];
   end1indexj=[]; end2indexj=[]; 
   [tdiff1, tdiff2, start1indexj, start2indexj, ...
    end1indexj, end2indexj, i, j] = tsub(i, j, shotcoord, fbcoord, fbtime);
    
   mddiff1=[]; deltf1=[]; delts1=[]; cvp1=[];
   [mddiff1, deltf1, delts1, cvp1] = pickTD(tdiff1, ...
       fbcoord(j,start1indexj:end1indexj), window, windmn, nd, ...
       offsetrange1, offsetrange2, shotcoord(i));
       
   mddiff2=[]; deltf2=[]; delts2=[]; cvp2=[];
   [mddiff2, deltf2, delts2, cvp2] = pickTD(tdiff2, ...
       fbcoord(j,start2indexj:end2indexj), window, windmn, nd, ...
       offsetrange1, offsetrange2, shotcoord(j));
   if(~isempty(cvp1))
      cvpi(i,j) = cvp1;
   end
   if(~isempty(cvp2))
      cvpj(i,j) = cvp2;
   end
end 
