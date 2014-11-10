function moveline(action)
% MOVELINE ... a simple utility to move a line around
%
% moveline
%
% MOVELINE is designed to allow the user to click on a line object in a
% figure window and drag it around. The mouse click can be with any button
% and only line objects can be moved. This is mainly a teaching tool to
% demo how to write interactive code.
%
% moveline ... after creating a figure window with at least one line object
%               in it, just type moveline and then click and drag your 
%               line(s) around.
%
% G.F. Margrave, CREWES, 2010
%

if(nargin<1) action='init'; end
if(strcmp(action,'init'))
  set(gcf,'windowbuttondownfcn','moveline(''click'')');
  return;
end
if(strcmp(action,'click'))
    hline=gco;
    obj=get(gco,'type');
    if(~strcmp(obj,'line'))
       msgbox('Click on a line!!!')
       return;
    end

    pt=get(gca,'currentpoint');
    set(hline,'userdata',pt(1,1:2));
    %set(hline,'erasemode','xor','linestyle','.');

    set(gcf,'windowbuttonmotionfcn','moveline(''move'')');
    set(gcf,'windowbuttonupfcn','moveline(''fini'')');
    return;
end
if(strcmp(action,'move'))
    hline=gco;

    pt1=get(hline,'userdata');
    pt2=get(gca,'currentpoint');
    pt2=pt2(1,1:2);

    del=pt2-pt1;

    x=get(hline,'xdata');
    y=get(hline,'ydata');

    set(hline,'xdata',x+del(1));
    set(hline,'ydata',y+del(2));
    set(hline,'userdata',pt2);

    return;
end

if(strcmp(action,'fini'))
    hline=gco;
    x=get(hline,'xdata');
    y=get(hline,'ydata');
    xmax=max(x);xmin=min(x);
    ymax=max(y);ymin=min(y);
    set(hline,'erasemode','normal','linestyle','-');
    set(hline,'marker','none');
    set(gcf,'windowbuttondownfcn','moveline(''click'')');
    set(gcf,'windowbuttonmotionfcn','');
    set(gcf,'windowbuttonupfcn','');
    %set(gca,'xlim',xlims);
    %set(gca,'ylim',ylims);
    xlims=get(gca,'xlim');
    ylims=get(gca,'ylim');
    if(xmin<xlims(1) || xmax>xlims(2))
        xint=round(diff(xlims)/100)*10;
        xmin=floor(xmin/xint)*xint;
        xmax=ceil(xmax/xint)*xint;
        set(gca,'xlim',[min([xlims(1) xmin]) max([xlims(2) xmax])]);
    end
    if(ymin<ylims(1) || ymax>ylims(2))
        yint=round(diff(ylims)/100)*10;
        ymin=floor(ymin/yint)*yint;
        ymax=ceil(ymax/yint)*yint;
        set(gca,'ylim',[min([ylims(1) ymin]) max([ylims(2) ymax])]);
    end

    return;
end	
