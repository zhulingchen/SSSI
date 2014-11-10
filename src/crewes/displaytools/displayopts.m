function displayopts(hfig)
% DISPLAYOPTS is a function that adds a contextmenu of display options to
% the figure.
%displayopts(hfig);
%
%hfig is the figure that you would like to add the display options menu to
% the options that are available in this function are
%               Increase Font Size
%               Decrease Font Size
%               Bold Font
%               Normal Font
%               Increase Line Weight
%               Decrease Line Weight
%               White Background
%               Grey Background
%               Large Figure (see bigfig for more information)
%               Small Figure (resizes the figure to the size that is was 
%                                when displayopts was called)
%               All Display Options ON              
%               All Display Options OFF
%               Copy Figure to Clipboard
%
% if no figure is supplied hfig defaults to gcf
%
% H.J.E Lloyd
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

if nargin<1
    hfig=gcf;
end
uni=get(hfig,'Units');
set(hfig,'Units','Normalized')
pos=get(hfig,'Position');
set(hfig,'Units',uni);

%checks to see if there is already a contextmenu in the figure then adds a
%new menu to the contextmenu or creates one
conmenu=get(hfig,'UIContextMenu');

if isempty(conmenu)
    conmenu=uicontextmenu;
    menus=conmenu;
else
    menus=uimenu('parent',conmenu, 'Label','Display Options','Separator','on');
end

%creates the menu options 
    uimenu('parent',menus, 'Label','Increase Font Size','Separator', 'off',...
        'Callback','bigfont(gcf,1.125,1)');
    uimenu('parent',menus, 'Label','Decrease Font Size','Separator', 'off',...
        'Callback','bigfont(gcf,(1/1.125),1)');
    uimenu('parent',menus, 'Label','Bold Font','Separator', 'off',...
        'Callback','bigfont(gcf,1,2)');
    uimenu('parent',menus, 'Label','Normal Font','Separator', 'off',...
        'Callback','bigfont(gcf,1,1)');
    uimenu('parent',menus, 'Label','Increase Line Weight','Separator', 'on',...
        'Callback','boldlines(gcf,2,1)');
    uimenu('parent',menus, 'Label','Decrease Line Weight','Separator', 'off',...
        'Callback','boldlines(gcf,(1/2),1)');
    uimenu('parent',menus, 'Label','White Background','Separator', 'on',...
        'Callback','whitefig');
    uimenu('parent',menus, 'Label','Grey Background','Separator', 'off',...
        'Callback','greyfig');
    uimenu('parent',menus, 'Label','Large Figure','Separator', 'off',...
        'Callback','bigfig(gcf)');
    uimenu('parent',menus, 'Label','Small Figure','Separator', 'off',...
        'Callback','smallfig(gcf,get(gcbo,''Userdata''))','Userdata',pos);
%     uimenu('parent',menus, 'Label','Zoom','Separator', 'on',...
%         'Callback','zoom');
%     uimenu('parent',menus, 'Label','Zoom Out','Separator', 'off',...
%         'Callback','zoom out;zoom off'); 
    uimenu('parent',menus, 'Label','All Display Options ON','Separator', 'on',...
        'Callback','bigfig(gcf);boldlines(gcf,4,2);bigfont(gcf,2,2);whitefig');
    uimenu('parent',menus, 'Label','All Display Options OFF','Separator', 'off',...
        'Callback',['smallfig(gcf,get(gcbo,''Userdata''));boldlines(gcf,.25,.5);',...
        'bigfont(gcf,.5,1);greyfig'],'Userdata',pos);
    uimenu('parent',menus, 'Label','Copy Figure to Clipboard','Separator', 'on',...
        'Callback','print(gcf,''-dbitmap'')');
    
    % adds the modified contextmenu back into the figure
    set(hfig,'uicontextmenu',conmenu)

end

