function [zcs,tcs]=readcheckshot(filename)
% [zcs,tcs]=readcheckshot(filename)
%
% Read in checkshot information from a standar ascii file. The file is
% defined as one with an arbitrary number of header lines containing any
% information followed by two columns of numbers: the first being depth
% and the second time. The header lines must begin with two slashes (//)
% but are otherwise arbitrary and are ignored by this function. The times
% may be one-way or two-way since CSCORRDISP is designed to handle both.
% Time will be output in seconds. If the time in the file contain any values
% greater than 10.0, they are assumed to be milliseconds and converted.
%
% G.F. Margrave November 1994
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
fid=fopen(filename);
if(fid==-1)
	zcs=-1;
	tcs=[];
	return;
end
%read past the header
p1=0;
buf=fgetl(fid);
p2=ftell(fid);
while(strcmp(buf(1:2),'//'))
	p1=p2;
	buf=fgetl(fid);
	p2=ftell(fid);
end
%reposition the file pointer
fseek(fid,p1,-1);
%read the data
dat=fscanf(fid,'%g');
zcs=dat(1:2:length(dat));
tcs=dat(2:2:length(dat));
zcs=zcs(:);
tcs=tcs(:);
tmax=max(tcs);
if(tmax>10)
	tcs=tcs/1000.;
end
fclose(fid);
