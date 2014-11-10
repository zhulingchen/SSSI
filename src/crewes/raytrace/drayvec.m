function drdt=drayvec(t,r)
% DRAYVEC: compute the derivative of ray vector (for vxz raytracing)
%
% drdt=drayvec(t,r)
%
% This is the fundamental driver for v(x,z) raytracing.
% The velocity model is defined by first calling rayvelmod.
% t ... scalar giveg the current time
% r ... 4x1 column vector giving (x,z,p,q)
% drdt ... 4x1 column vector giving the time-derivative of r
%
% By default, DRAYVEC uses nearest neighbor interpolation. Bilinear 
% interpolation is available for more accuracy. To get bilinear 
% interpolation, define a global variable called BILINEAR (all caps) 
% and set its value to 1. This is quite a bit slower than nearest 
% neighbor for the same grid size.
%
%
% by G.F. Margrave, CREWES, June 2000
%
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.
 
% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by 
% its author (identified above) and the CREWES Project.  The CREWES 
% project may be contacted via email at:  crewesinfo@crewes.org
% 
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) Use of this SOFTWARE by any for-profit commercial organization is
%    expressly forbidden unless said organization is a CREWES Project
%    Sponsor.
%
% 2) A CREWES Project sponsor may use this SOFTWARE under the terms of the 
%    CREWES Project Sponsorship agreement.
%
% 3) A student or employee of a non-profit educational institution may 
%    use this SOFTWARE subject to the following terms and conditions:
%    - this SOFTWARE is for teaching or research purposes only.
%    - this SOFTWARE may be distributed to other students or researchers 
%      provided that these license terms are included.
%    - reselling the SOFTWARE, or including it or any portion of it, in any
%      software that will be resold is expressly forbidden.
%    - transfering the SOFTWARE in any form to a commercial firm or any 
%      other for-profit organization is expressly forbidden.
%
% END TERMS OF USE LICENSE

global RTV2 RTDLNVDX RTDLNVDZ RTDG BILINEAR RTX RTZ 
if(isempty(RTV2))
   error('velocity model not defined. Use RAYVELMOD')
end
if(isempty(BILINEAR))
   BILINEAR=0;
end

[m,n]=size(RTV2);

x=r(1);z=r(2);p=r(3);q=r(4);

drdt=zeros(length(r),1);

%determine current position
xx=(x-RTX(1))/RTDG+1; zz=(z-RTZ(1))/RTDG+1;% Actual fractional grid
% guard against leaving model
xx=min([xx n-1]);xx=max([xx 1]);
zz=min([zz m-1]);zz=max([zz 1]); 

if(~BILINEAR)
   %nearest neighbor interpolation
   v2=RTV2(round(zz),round(xx));
   dlnvdx=RTDLNVDX(round(zz),round(xx));
   dlnvdz=RTDLNVDZ(round(zz),round(xx));
else


   %bilinear interpolation of v squared

   ix=floor(xx);%Upper left
   iz=floor(zz);%Upper left

   factx=xx-ix;
   factz=zz-iz;
   factxz=factx*factz;
   a=RTV2(iz,ix) + factx*(RTV2(iz,ix+1)-RTV2(iz,ix));
   b=RTV2(iz+1,ix) + factx*(RTV2(iz+1,ix+1)-RTV2(iz+1,ix));
   v2= a+factz*(b-a);

   %bilinear interpolation of logderivs
   a=RTDLNVDX(iz,ix) + factx*(RTDLNVDX(iz,ix+1)-RTDLNVDX(iz,ix));
   b=RTDLNVDX(iz+1,ix) + factx*(RTDLNVDX(iz+1,ix+1)-RTDLNVDX(iz+1,ix));
   dlnvdx=a+factz*(b-a);
   a=RTDLNVDZ(iz,ix) + factx*(RTDLNVDZ(iz,ix+1)-RTDLNVDZ(iz,ix));
   b=RTDLNVDZ(iz+1,ix) + factx*(RTDLNVDZ(iz+1,ix+1)-RTDLNVDZ(iz+1,ix));
   dlnvdz=a+factz*(b-a);

end

drdt(1) = v2*p;
drdt(2) = v2*q;
drdt(3) = -dlnvdx;
drdt(4) = -dlnvdz;
