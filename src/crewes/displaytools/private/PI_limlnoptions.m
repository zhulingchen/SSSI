function PI_limlnoptions();
h=get(gcf,'userdata');
hlimbox=h(14);
limdat=get(hlimbox,'userdata');
limdat={limdat{2} limdat{1} limdat{4}};
stdat=get(gcbo,'userdata');
st={'linestyle' 'marker' 'color'};
ls={'-' '--' '-.' ':'};
mk={'*' 'x' 'o' '+' '.'};
col={'r' 'g' 'b' 'k' [1 0 1]};
ch={ls mk col};
ch=ch{stdat(1)};
set(limdat{2},st{stdat(1)},ch{stdat(2)});
set(limdat{1},st{stdat(1)},ch{stdat(2)});
if(stdat(1)==3)
    set(limdat{3},st{stdat(1)},ch{stdat(2)});
end
