function setmenustatus( menulist, dataname )
if( ~isempty(refdata('get', dataname)) )
   flag = 'on';
else
   flag = 'off';
end
for i=1:length(menulist)
   set(menulist(i), 'enable', flag);
end
