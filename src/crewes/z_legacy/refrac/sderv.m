function[delts1,delts2,avgmdxs1,avgmdxs2]=sderv(mddiff3,mddiff4,avgmdx1,...
avgmdx2,delt1,delt2)
clear delts1
%clear dts2
[m,n]=size(avgmdx1);
if( size(mddiff3) ~= 0 )
  delts1=delt1(1:n-2)-delt1(3:n);
  avgmdxs1=(avgmdx1(1:n-2)+avgmdx1(3:n))/2;
%  dts1 = clip(delts1,1);
else
  avgmdxs1 = [];
end
clear delts2
%clear dts2
[m,n]=size(avgmdx2);
if( size(mddiff4) ~= 0 )
  delts2=delt2(1:n-2)-delt2(3:n);
  avgmdxs2=(avgmdx2(1:n-2)+avgmdx2(3:n))/2;
%  dts2 = clip(delts2,1);
else 
  avgmdxs2 = [];
end
