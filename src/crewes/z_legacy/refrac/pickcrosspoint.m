
function curpt(action)
 if(nargin<1) %set the button down function
                set(gcf,'windowbuttondownfcn','curpt(''init'')');
                return;
        end
        if(strcmp(action,'init'))
                %hline=gco;
                %if(~strcmp(get(hline,'type'),'line'))
                %       return;
                %end
                pt=get(gca,'currentpoint');
                %set(hline,'userdata',pt(1,1:2));
                %set(hline,'erasemode','xor','linestyle','.');
                %set(gcf,'windowbuttonmotionfcn','curpt(''move'')');
                %set(gcf,'windowbuttonupfcn','curpt(''fini'')');
                line(pt(1,1),pt(1,2),'color','r','linestyle','*');
                return;
        end
        if(strcmp(action,'move'))
                hline=gco;
                pt1=get(hline,'userdata');
