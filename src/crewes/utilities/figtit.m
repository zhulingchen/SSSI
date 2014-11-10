function figtit(str,figno)

if(nargin<2)
figno=gcf;
end

set(figno,'name',str)
