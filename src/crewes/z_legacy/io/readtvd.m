function [zmd,azmi,incl,utmx,utmy,utmzone]=readtvd(filename)
% [zmd,azmi,incl,utmx,utmy]=readtvd(filename)
%
% Read a well bore deviation file survey. The file is
% defined as consisting of any number of headers followed
% by three columns of ascii numbers: measured_depth,
% azimuth, inclination. Headers are denoted by // in the
% first two columns.
%
% Note, the first header line is examined for well surface
% location information. If found, these are returned as
% utmx, utmy and utmzone
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
%open the file
fid=fopen(filename);
rec=fgetl(fid);
utmx=0;
utmy=0;
utmzone=0;
%read past comments
while( strcmp(rec(1:2),'//') )
 %look for location info
 ind=findstr(rec,'X=');
 if(~isempty(ind))
		ind2=findstr(rec,'Y=');
		ind3=findstr(rec,'ZONE=');
		utmx=sscanf(rec(ind+2:length(rec)),'%lg');
		utmy=sscanf(rec(ind2+2:length(rec)),'%lg');
		utmzone=sscanf(rec(ind3+5:length(rec)),'%d');
	end
	rec=fgetl(fid);
end
%load the rest of the file
dat=fscanf(fid,'%g',[3,inf]);
fclose(fid);
dat=dat.';
zmd=dat(:,1);
azmi=dat(:,2);
incl=dat(:,3);
