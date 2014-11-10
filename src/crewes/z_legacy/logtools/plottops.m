function htops=plottops(tops,topnames,klr)
% PLOTTOPS: plot tops on a log plot
%
% htops=plottops(tops,topnames,klr)
%
% PLOTTOPS will plot horizontal lines across the current figure
% for each of a set of formation tops. The tops will be labeled
% and will extend for the length of the x axis. It is assumed
% that the log is plotted with the y coordinate being depth
% (or time).
% 
% tops ... vector of vertical coordinates (depth or time as 
%	appropriate for the figure) for the tops
% topnames ... string matrix containing the top names. One name
%	per row.
% klr ... color to plot with
% ************* default = 'k' (black) ************
% htops ... vector of handles of the tops and their text labels.
% 	They can be deleted by the command: delete(htops).
%
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.
%
% G.F. Margrave, CREWES, 1995
%
 
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
if(nargin<3)
    klr='k';
end
if(~isstr(tops))
	action='plot';
else
	action=tops;
end
if(strcmp(action,'plot'))
	%determine the x limits
	xlim=get(gca,'xlim');
	%col = [1 .5 0];
    
%     %find a log  in the axis
%     hkids=get(gca,'children');
%     hlog=0;
%     for kkk=1:length(hkids)
%         if(strcmp(get(hkids(kkk),'type'),'line'))
%             hlog=hkids(kkk);
%         end
%     end
%     if(hlog==0)
%         error('plot a log first then plot tops')
%     end
%     
%  xlog=get(hlog,'xdata');
%  ylog=get(hlog,'ydata');
%  c=polyfit(ylog,xlog,3);
%  polylog=polyval(c,ylog);
 ntops=length(tops);
	htops=zeros(3*ntops+1,1);
	for k=1:ntops
		if( tops(k) > 10 )
			lbl=sprintf('%5.2f',tops(k));
		else
			lbl=sprintf('%2.5f',tops(k));
		end
        %indy=near(ylog,tops(k));
		%htops(ntops+k)=text(min([polylog(indy(1))+.25*(xlim(2)-xlim(1)) xlim(2)]),tops(k),topnames(k,:),...
        htops(ntops+k)=text(min(xlim(2)),tops(k),topnames(k,:),...
			'verticalalignment','baseline',...
			'horizontalalignment','left',...
			'color',klr,'fontsize',9);
        
        %htops(2*ntops+k)=text(max([polylog(indy(1))-.25*(xlim(2)-xlim(1)) xlim(1)]),tops(k),topnames(k,:),...
        htops(2*ntops+k)=text(xlim(1),tops(k),topnames(k,:),...
			'verticalalignment','baseline',...
			'horizontalalignment','right',...
			'color',klr,'fontsize',9);
		htops(k)=line(xlim,[tops(k) tops(k)],...
			'color',klr,'linestyle','-.',...
			'userdata',['Top: ' topnames(k,:) ' value: ' lbl],...
			'buttondownfcn','plottops(''sayhey'')');
	end
	hmsg=uicontrol('style','text','units','normalized','position',[0 0 1 .05],...
		'string','Click MB3 on any top to see its name and value');
	set(gcf,'userdata',hmsg);
	htops(2*ntops+1)=hmsg;
	return;
end
if(strcmp(action,'sayhey'))
	hmsg=get(gcf,'userdata');
	flag= get(gcf,'selectiontype');
	if(~strcmp(flag,'alt'))
		set(hmsg,'string',...
			'Click MB3 on any top to see its name and value');
		return;
	end
	msg=get(gco,'userdata');
	set(hmsg,'string',msg);
	return;
end
