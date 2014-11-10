function greyfig
% GREYFIG the figure background to grey
%
% GREYFIG changes the current figure's default background color
%          to grey.
%
%	G.F. Margrave
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

set(gcf,'color',[.8 .8 .8])
set(gcf,'defaulttextcolor','black')
set(gcf,'inverthardcopy','on')

%find all axes
% hkids=get(gcf,'children');
% axid=zeros(size(hkids));
% for k=1:length(hkids)
% 	if(strcmp(get(hkids(k),'type'),'axes'))
% 		axid(k)=1;
% 	end
% end
% 
% haxes=hkids(find(axid==1));
% for kax=1:length(haxes)
% 	hax=haxes(kax);
% 	set(hax,'xcolor','black')
% 	set(hax,'ycolor','black')
% 	set(hax,'zcolor','black')
% 	htit=get(hax,'title');
% 	if( htit>0 )
% 		set(htit,'color','black');
% 	end
% 	h=get(hax,'children');
% 	for k=1:length(h)
% 		if( strcmp(get(h(k),'type'),'text') |  ...
% 			strcmp(get(h(k),'type'),'Text') )
% 			kol=get(h(k),'color');
% 			if( sum(kol)==3 )
% 				set(h(k),'color','black');
% 			end
% 		end
% 		if( strcmp(get(h(k),'type'),'line') | ...
% 			strcmp(get(h(k),'type'),'Line') )
% 			kol=get(h(k),'color');
% 			if( sum(kol)==3 )
% 				set(h(k),'color','black');
% 			end
% 		end
% 	
% 	end
% 
% end
% 
