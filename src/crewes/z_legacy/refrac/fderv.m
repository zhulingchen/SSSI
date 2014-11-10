function[delt1,delt2,avgmdx1,avgmdx2]=fderv(mddiff3,mddiff4,x1,x2)
clear delt1
%clear dt1
[m,n]=size(x1);
if( size(mddiff3) ~= 0 )
  delt1=mddiff3(1:n-2)-mddiff3(3:n);
  avgmdx1=(x1(1:n-2)+x1(3:n))/2;
%  dt1 = clip(delt1,1);
else
  avgmdx1 = [];
end
clear delt2
%clear dt2
[m,n]=size(x2);
if( size(mddiff4) ~= 0 )
  delt2=mddiff4(1:n-2)-mddiff4(3:n);
  avgmdx2=(x2(1:n-2)+x2(3:n))/2;
%  dt2 = clip(delt2,1);
else 
  avgmdx2 = [];
end
