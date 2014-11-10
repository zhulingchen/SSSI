function lineinfo=newlinefini
%
% lineinfo=newlinefini
% 
% NEWLINEFINI is called to obtain the handle and (x,y)'s of a line created
% by NEWLINE and to perform needed cleanup. The returned information
% consists of the handle of the new line, a vector of x coordinates and a
% vector of y coordinates. Note that the number of points on the new line is
% not provided but is simply: n=(length(lineinfo)-1)/2
%
% lineinfo=[key handle x1 x2 x3 .... xn y1 y2 y3 .... yn]
%
% The key gives needed information about the new line as follows:
% key = 1 ... the line was created by drawing with the puck
% key = 2 ... the line was created by duplicating an existing line
% key = 3 ... the line was created by undoing a line deletion
% key = 4 ... the line was actually deleted. Its handle is no longer valid.
%
% When a line is deleted, the lineinfo return value has all of the above
% information plus additionaly info encoding the deleted lines color,
% linestyle, etc. The is appended to the above information using a NaN as
% a separator. Thus if key==4, then :
%   ind=find(isnan(lineinfo));
%	info = lineinfo(1:ind-1);
% will store in info the simple lineinformation vector given above. The
% simplest way to restore a deleted line with its color etc intact is to call:
% newline('undo',lineinfo,transfer);
% where lineinfo is the complete return value from NEWLINEFINI and transfer
% is a string indicatining where the NEWLINE should return to when completed
% (see NEWLINEINIT)
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
hstor=[]; hundo=[];
% clean up graphics if needed
	% Get the storage buckets
	h=get(gcf,'children');
	found=0;
	for k=1:length(h)
		if( strcmp(get(h(k),'type'),'uicontrol') )
			if( strcmp(get(h(k),'style'),'text') )
				if( strcmp(get(h(k),'string'),'yrag_new') )
					hstor=h(k);
					found=found+1;
				end
				if( strcmp(get(h(k),'string'),'yrag_trans') )
					hundo=h(k);
					found=found+1;
				end
				if( found== 2)
					break;
				end
			end
		end
	end
% test userdata of hstor. if it is null, then a normal exit was achieved, else, it
% contains the information we want to return
if(~isempty(hstor))
	stuff=get(hstor,'userdata');
	if(~isempty(stuff))
		lineinfo=stuff;
	else
	 lineinfo=get(gca,'userdata');
	end
else
	lineinfo=[];
end
	
% delete the storage buckets
	delete(hstor);
	delete(hundo);
	
	set(gca,'userdata','');
	%set the buttondown function to null
	set(gcf,'windowbuttondownfcn','');
