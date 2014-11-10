function code=writeseistvd(filename,x,y,z,zmd,utmx,utmy,uwid,tz,zt,flag)
% code=writeseistvd(filename,x,y,z,zmd,utmx,utmy,uwid,tz,zt,flag)
%
% WRITESEISTVD writes a SEISLINE readable ascii file containing
% the well deviation geometry in UTM coordinates
%
% filename ... full filename of the file to write to
% x,y,z ... 3 column vectors giving the borehole geometry
% zmd ... vector of same length as z giving the measured depths
% utmx ... utm x coordinate of well at z==0
% utmy ... utm y coordinate of well at z==0
% uwid ... unique well id (ugly)
% tz,zt ... time-depth function for well
% flag ... 1 -> zt (depths of time-depth function) are measured depths
%	along the borehole and must be converted to tvd's
%      ... 0 -> zt is already in tvd
%  ***** default for flag is 1 ****
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
if (nargin < 11)
	flag=1;
end
%open the file
fid=fopen(filename,'w');
if(fid==-1)
	code=-1;
	return;
end
% compute times to go with the depths
if(flag)
	t=interpextrap(zt,tz,zmd);
else
	t=interpextrap(zt,tz,z);
end
% compute the pretty well id
pid=prettywid(uwid);
%write first line
rec='CH: WELL NAME       SM    KB FLG CL DIST AZM SEIS_CUBE';
rec=[rec '  CHEVNO STK  CMP#   SURF_X   SURF_Y     TD'];
fprintf(fid,'%s\n',rec);
%second line
fmt='WH: %s  0     0   1  1  500   0    ***      ******   *    ** ';
fmt=[fmt '%8.0f %8.0f   %5.0f'];
rec=sprintf(fmt,pid,utmx,utmy,z(length(z)));
fprintf(fid,'%s\n',rec);
%third line
rec='CP: WELL NAME              WX       WY   XOFF   YOFF XLAB';
rec=[rec ' YLAB   TVD  TIME  WELL_ID'];
fprintf(fid,'%s\n',rec);
%now the rest
fmt='WP: %s   %8.0f %8.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f';
zip=0.;
uno=1.;
for k=1:length(z)
 ux=utmx+x(k);
 uy=utmy+y(k);
	rec=sprintf(fmt,pid,ux,uy,x(k),y(k),zip,zip,z(k),1000*t(k),uno);
	fprintf(fid,'%s\n',rec);
end
	
fclose(fid)
