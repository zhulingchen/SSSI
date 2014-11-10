function [mnems,units,values,descs,cmnts,irc]=lashblkread(block,verflag)
% [mnems,units,values,descs]=lashblkread(block,verflag)
%
% LASHBLKREAD interrogates an LAS logical data block (see lashgetblk)
% to decompose it into more elemental string matricies. If one of the
% string matricies is known to actually be numerical, it can be
% simply converted with str2num.
% block ... string matrix containing the parameter block
% verflag ... LAS version number. Defaults to 2.0
% mnems ... string marix of the first field on each line. Usually this is
%	a four character mnemonic. In the case of a ~tops block, it is a more
%	lengthy formation name
% units ... string matrix of physical units. May be blank or null
% values ... string matrix of values. These often hav a numerical
%	interpretation.
% descs ... string matrix of colloquial descriptions
% cmnts ... string matrix of comments
% irc ... vector of row numbers on which the comments appear in the
%	block.
%
% G.F. Margrave, Department of Geology and Geophysics,
%	University of Calgary, 1996
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
if (nargin<2) verflag=2; end
%special code for tops
if(strcmp(lower(block(1,1:5)),'~tops'))
	k=2; %skip line 1
	[nlines,m]=size(block);
	mnems=[];
	units=[];
	values=[];
	descs=[];
	cmnts=[];
	irc=[];
	while(k<=nlines)
		str=block(k,:);
		ls=length(str);
		indn32=find(abs(str~=32));%non blanks
		ind32=find(abs(str==32));%blanks
		if(~isempty(indn32))% test for blank lines and ignore
			i1=indn32(1);
		else
			i1=0;
		end
		if(i1)
			%test for comment
			if(str(i1)=='#')
				if(isempty(cmnts))
					cmnts=str;
				else
					cmnts=str2mat(cmnts,str);
				end
				irc=[irc k];
			else
				%search for .
				ind=findstr(str,' . ');
				if(isempty(ind))
					mnem=sscanf(str,'%s');
					val=sscanf(str,'%*s %s');
				else
					iv=find(indn32<ind);
					iv=indn32(iv);
					mnem=str(min(iv):max(iv));
					iv=find(indn32>ind+2);
					iv=indn32(iv);
					val=str(min(iv):max(iv));
					%test for a :
					ind=find(val==':');
					if(~isempty(ind))
						val=val(1:ind(1)-1);
					end
				end
				if(isempty(mnems))
					mnems=mnem;
					values=val;
				else
					mnems=str2mat(mnems,mnem);
					values=str2mat(values,val);
				end
				
			end
		end
		
		k=k+1;
	end
	[n,m]=size(mnems);
	units=setstr(32*ones(n,1));
	descs=setstr(32*ones(n,1));
	return;
end
if(verflag>=1.2)
	k=2; %skip first line because its just the block title
	[nlines,m]=size(block);
	mnems=[];
	units=[];
	values=[];
	descs=[];
	cmnts=[];
	irc=[];
	while(k<=nlines)
		str=block(k,:);
		ls=length(str);
		indn32=find(abs(str~=32));%non blanks
		ind32=find(abs(str==32));%blanks
		if(~isempty(indn32))% test for blank lines and ignore
			i1=indn32(1);
		else
			i1=0;
		end
		if(i1)
			%test for comment
			if(str(i1)=='#')
				if(isempty(cmnts))
					cmnts=str;
				else
					cmnts=str2mat(cmnts,str);
				end
				irc=[irc k];
			else
				%find the mnem
				ind=find(str=='.');
				if(isempty(mnems))
					mnems=str(i1:ind(1)-1);
				else
					mnems=str2mat(mnems,str(i1:ind(1)-1));
				end
				
				%find the units
				in1=between(ind(1),10,indn32);
				if(~in1)
					un=' ';
					iv1=11;%start looking for values here
				else
					i1=indn32(in1(1));%first nonblank > '.' and < 10
					in2=find(ind32>i1);
					i2=ind32(in2(1))-1;
					un=str(i1:i2);
					
					iv1=i2+1;
				end
				if(isempty(units))
					units=un;
				else
					units=str2mat(units,un);
				end
				
				%find the values
				icolon=find(str(iv1:ls)==':');
				if(isempty(icolon))
					icolon=ls+1;
				else
					icolon=iv1+icolon-1;
				end
				in1=between(iv1(1),icolon(1)-1,indn32);
				if(in1==0)
					val=' ';
				else
					i1=indn32(in1(1));
					i2=indn32(in1(length(in1)));
					val=str(i1:i2);
				end
				if(isempty(values))
					values=val;
				else
					values=str2mat(values,val);
				end
				
				%the decsriptions
				if(icolon>ls)
					des=' ';
				else
					des=str(icolon+1:ls);
				end
				if(isempty(descs))
					descs=des;
				else
					descs=str2mat(descs,des);
				end
				
				
			end
		end
		
		k=k+1;
	end
else
	error('lasblkread can''t handle earlier than version 1.2 files')
end
		
	
