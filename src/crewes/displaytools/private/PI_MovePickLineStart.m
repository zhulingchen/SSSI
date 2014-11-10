function PI_MovePickLineStart();
nm=get(gco,'tag');
switch nm
    case 'PICKMARKER'
        set(gcf,'name','Click & Hold MB1 on Markers to move.  MB3 menu, or function change to stop line moving')
        set(gcf,'windowbuttonmotionfcn','plotimage(''MovePickLineMotion'')',...
            'windowbuttonupfcn','plotimage(''MovePickLineEnd'')');
        udat=get(gco,'userdata');
        set(gco,'erasemode','xor');
        set(findobj(gcf,'type','line','tag','PICKS'),'erasemode','xor');
    case 'PICKS'
        sltype=get(gcf,'selectiontype');
        switch sltype
            case 'alt'
                % will continue on and open line menu
            case 'normal'
                udat=get(gco,'userdata');
                if(isempty(udat))
                elseif(~ishandle(udat(1)))
                else
                    cpt=get(gca,'currentpoint');
                    set(gca,'userdata',cpt);
                    set(gco,'erasemode','xor');
                    set(udat(1),'erasemode','xor');
                    set(udat(2),'erasemode','xor');
                    set(gcf,'name','Click & Hold MB1 on Markers to move.  MB3 menu, or function change to stop line moving')
                    set(gcf,'windowbuttonmotionfcn','plotimage(''MovePickLineMotion'')',...
                        'windowbuttonupfcn','plotimage(''MovePickLineEnd'')');
                    set(findobj(gcf,'type','line','tag','PICKS'),'erasemode','xor');
                end
            case 'More To Come'
        end    
end
