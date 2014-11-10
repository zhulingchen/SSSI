function h=clickaline(action)
% CLICKALINE ... create a line in the current axis by mouse clicks
%
% clickaline
%
% To use this function, just issue the command: clickaline
%
% Usage: MB1: press mouse button 1 to define a point on a line
%        MB3: mouse button 3 signals done (and does not add a point)
%        MB2: mouse button 2 on any object of type 'line' deletes it
% When completed (as signaled by MB3) the user data of the current axis
% contains a structure called linespec with the fields:
% linespec.handle ... the handle of the created line
% linespec.xdata ... the x data of the created line
% linespec.ydata ... the y data of the created line
%
% by G.F. Margrave, CREWES, 2010
%
if(nargin<1)
    action='init';
end
kol='k';
switch action
    case 'init'
        set(gcf,'windowbuttondownfcn','clickaline(''newpt'')');
        set(gca,'userdata',[]);
    case 'newpt'
        pt=get(gca,'currentpoint');
        linespec=get(gca,'userdata');
        button=get(gcf,'selectiontype');
        if(strcmp(button,'normal'))
            if(isempty(linespec))
                h=line(pt(1,1),pt(1,2),'marker','*','color',kol);
                linespec.handle=h;
                linespec.xdata=pt(1,1);
                linespec.ydata=pt(1,2);
                set(gca,'userdata',linespec);
            else
                xdata=[linespec.xdata pt(1,1)];
                ydata=[linespec.ydata pt(1,2)];
                set(linespec.handle,'xdata',xdata,'ydata',ydata);
                linespec.xdata=xdata;
                linespec.ydata=ydata;
                set(gca,'userdata',linespec);
            end
        elseif(strcmp(button,'alt'))
            set(gcf,'windowbuttondownfcn','');
        elseif(strcmp(button,'extend'))
            if(strcmp(get(gco,'type'),'line'))
                delete(gco);
            end
        end
end
            