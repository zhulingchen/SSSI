function [X,Y,NAME] = readwellfile( filename )
% [X,Y,NAME] = readwellfile( filename )
%
% where X and Y hold coordinates (UTM) of the well positions 
% where NAME holds the names of the wells (12 char). 
% where filename is the ascii card file , where each line
%   holds the x,y,name of each well.  (This file is made using
%   readwells, a program made in the seisline environment).
%
% by T.N. Bishop, November 1993
%   see also seis2well, linewelltie
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
%filename = 'wells.dat';
[fid,message] = fopen(filename,'r');
if( fid < 0 )
	disp(' Unable to open file because:');
	disp( message );
	error(' file I/O error in readwellfile');
end
% read the next line and get the character string s 
s=fgetl(fid);
[message,errnum] = ferror(fid);
X = [];
Y = [];
NAME = [];
while   errnum == 0
  X = [X;sscanf(s,'%f',1)];
  Y = [Y;sscanf(s,'%*f %f',1)];
  nn = [sscanf(s,'%*f %*f %9s',1) zeros(1,12)];
  NAME = [NAME;nn(1:12)];
  s=fgetl(fid);
  [message,errnum] = ferror(fid);
end
message;
errnum;
disp(['number of wells = ',num2str(length(Y))]);
