function makehor(t, utmx, utmy, shotpt, hor, linename)
   
% makehor(t, utmx, utmy, shotpt, hor, linename)
%   supplements makeline2
%
%   Tom Bishop, CPTC
%   December 8, 1994
%
% Note: makehor is used by LOGSEC to write out a horizons file which can be
% taken into SEISLINE. The file name is hardwired to be: logsec_hor.dat
% Once this file is generated it is run through John Bornhurst's program
% er2slss.csh which generates 2 files for SEISLINE: xc.hor and 
% survey_summary.hor . Meanings of the variables are:
%	t ... matrix of times for horizons one per seismic trace.
%		One horizon per row. Times in milliseconds,
%		use -9999. to flag missing values.
%	utmx, utmy ... utm coordinates for the traces.
%	shotpt ... shotpoints
%	linename ... string containg the linename. Should be the same as
%		used in makeline2
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

% warning, this will write over any existing file named
% logsec_hor.dat
filename='logsec_hor.dat';
[fid,message] = fopen(filename,'w');
if( fid < 0 )
    disp('problem opening file named logsec.dat');
	disp( message );
	error(' file I/O error ');
end

%print header
[nhor,ntrace]=size(t);
for i = 1:nhor
	fprintf(fid,'! hor=%20s, att%2d = time\n',hor(i,:),i);
end

%print trace raster picks, extract_raster format
for i = 1:ntrace
    for j = 1:nhor
		if( t(j,i) > -9999. & t(j,i) < 999999. )
		  AT(j,1:7) = sprintf('%7.1f',t(j,i));
		else
	      AT(j,1:7) = '       ';
		end
	end
	XX(1:10) = sprintf('%10.2f',utmx(i));
	YY(1:11) = sprintf('%11.2f',utmy(i));
	SP(1:8) = sprintf('%8.2f',shotpt(i));
    fprintf(fid,'%-12s%-8s%-10s%-11s', linename,SP(:),XX(:),YY(:) );
    for j = 1:nhor
      fprintf(fid,'%7s',AT(j,:));
	end
    fprintf(fid,'\n');
end
fclose(fid);
disp(' ')
disp(' File named logsec_hor.dat was successfully created.');
disp(' It is an ascii file in extract_raster zmap format.'); 
disp(' It can be read by the program er to be converted');
disp('    to SEISLINE .hor files.');

%  below is format used in extract raster for reference:
%       write(2,330)event(jj),natt
%330    format('! hor=',A20,', att',I2,' = time')
%		write(2,'(a12,a8,F10.2,F11.2,6A7,I6)')
%		*         line_name,sp_name,X,Y,AT(1),AT(2),AT(3),
%		*         AT(4),AT(5),AT(6),index
