function writelas(fullfilename,lasheader,logmat,z)

% writelas(fullfilename,lasheader,logmat,z)
% writelas(fullfilename,lasheader,logmat)
%
% WRITELAS writes and LAS (Log ASCII Standard) dataset to disk.
%
% fullfilename ... a string containing the fully qualified filename
% lasheader ... a string matrix containing as las header which
%		matches the data in z and logmat. Two easy ways to get this
%		are 1) from READLAS or 2) from MAKELASHEADER
% logmat ... matrix of logs, one per column. Note that since LAS
%	does not support NaN's, all NaNs in the logmatrix must be set
%	to a NULL value consistent with that in the LAS header.
% z ... vector of depths
%
% NOTE: if Z is not supplied then it is assumed that the first column
%  of logmat containes the depths. If z is supplied, the it is
%  prepended to logmat as column 1.
%
% NOTE: number of rows in logmat must be the same as the length
% of z
%
% G.F. Margrave
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

if( nargin==4 )
	logmat=[z(:) logmat];
end

%test for computer. Special method for mac
comp=computer;
if(strcmp(comp,'MAC2')|strncmp(comp,'PCWIN',5))
	fid=fopen(fullfilename,'w');
	if(fid==-1)
	   msgbox('Cannot open output file. Perhaps it is marked as read only?')
	   return
	end
	[nlines,nchar]=size(lasheader);
	for kline=1:nlines
		fprintf(fid,'%s\r\n',lasheader(kline,:));
	end
	%now build format string for the logmat
	[ns,nlogs]=size(logmat);
	fmtstr='';
	for k=1:nlogs
		fmtstr=[fmtstr '%13.7e   '];
	end
	fmtstr=[fmtstr '\r\n'];
	
	%write out the logmat
	fprintf(fid,fmtstr,logmat');
	
	fclose(fid);
else
	fid=fopen(fullfilename,'w');
	[nlines,nchar]=size(lasheader);
	for kline=1:nlines
		fprintf(fid,'%s\r\n',lasheader(kline,:));
	end
	fclose(fid);
	tmpname=['/tmp/junk' num2str(fix(clock)) '.tmp'];
	indies=find(tmpname==' ');
	tmpname(indies)='';
	eval(['save ' tmpname ' logmat -ascii']);
	unixcmd = ['cat ' tmpname ' >> ' fullfilename];
	unix(unixcmd);
end

disp(['file ' fullfilename ' made ']);
