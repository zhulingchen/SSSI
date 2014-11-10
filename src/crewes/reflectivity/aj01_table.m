%ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ            
%                                                                               
%                  FUNCTION AJ0                                                 
%                                                                               
%ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ          
%                                                                               
function [aj0, aj1] = aj01_table(z,zmax,dz)
%
%ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
%                                                                               
%   AJ01 CALCULATES THE VALUE  OF A BESSEL FUNCTION OF THE FIRST KIND OF         
%   ORDERS 0 AND 1 FOR REAL ARGUEMENTS.        
%                                                                               
      a=[0.0002100,-0.0039444,0.0444479,-0.3163866,...                         
            1.2656208,-2.2499997];                                              
      b=[0.00001109,-0.00031761,0.00443319,-0.03954289,...                     
            0.21093573,-0.56249985];                                          
%
%ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
%
persistent besseltable zmaxtable dztable
if(isempty(besseltable))
    maketable=1;
elseif(zmaxtable~=zmax||dztable~=dz)
    maketable=1;
else
    maketable=0;
end
if(maketable==1)
      zz=0:dz:zmax;
      zmaxtable=zmax;
      dztable=dz;
      besseltable=zeros(2,length(zz));
      time1=clock;
      disp('making Bessel table');
      for k=1:length(zz)
          z=zz(k);
          if ( abs(z) < 3.0) 
             aj0 = 0.0;
             aj1 = 0.0; 
             x2 = (z/3.0)^2;                                                            
             for i = 1 : 6                                                                
                aj0 = ( aj0 + a(i) )*x2;
                aj1 = ( aj1 + b(i) )*x2;
             end                                                                
             aj0 = aj0 + 1.0;
             aj1 = ( aj1 + 0.5 )*z;
          else                                                                   
             sqrx = sqrt(2.0 / (pi * z));
             aj0 = sqrx * cos(z - pi * 0.25);
             aj1 = sqrx * cos(z - pi * 0.75);
          end  
          besseltable(1,k)=aj0;
          besseltable(2,k)=aj1;
      end
      time2=clock;
      timeused=etime(time2,time1);
      disp(['Bessel table completed in ' int2str(timeused) ' seconds'])
else
    nz=round(z/dztable)+1;
    aj0=besseltable(1,nz);
    aj1=besseltable(2,nz);
end

