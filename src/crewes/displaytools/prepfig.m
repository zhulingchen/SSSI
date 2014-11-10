bigfig; %enlarge the figure to get more pixels
bigfont; %enlarge the fonts in the figure
boldlines(gcf,4,2); %make lines and symbols "fatter"
whitefig; %make the background white
hideui; %hide any user interface controls
if(ispc)
    print -dbitmap %copy the figure to the clipboard (windows only)
else
    print -dtiff
end