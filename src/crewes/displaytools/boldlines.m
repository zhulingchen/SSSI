function boldlines(hdaddy,linewidth,markersize,kol)
% BOLDLINES: changes the thickness of lines and the saize of markers
%
% boldlines(hdaddy,linewidth,markersize)
%
% Its simple, just type "boldlines" at the matlab prompt. The lines are
% made much thicker so that the figures can make good slides.
%
% hdaddy ... handle of a figure or axes object or other parent object 
%		which contains line objects
%   ******** default = gcf ******
% linewidth ... desired line width expressed as a ratio of output
%		to input size
%   ******** default 4 ********
% markersize ... desired line width expressed as a ratio of output
%		to input size
%   ******** default 2 ********
% kol ... color of marker or line
%   ******** default is no color change ********
%
% You can customize the default behavior of boldlines by defining some
%  globals: BOLDLINES_LW and BOLDLINES_MS . Set these to have your desired
%  default values of linewidth and markersize. A good place to define
%  these is in your startup.m file.
%
% Gary Margrave, CREWES
%


if(nargin<3)
    global BOLDLINES_MS
    if(isempty(BOLDLINES_MS))
	    markersize=2;
    else
        markersize=BOLDLINES_MS;
    end
end
if(nargin<2)
    global BOLDLINES_LW
    if(isempty(BOLDLINES_LW))
	    linewidth=4;
    else
        linewidth=BOLDLINES_LW;
    end
end
if(nargin<1)
	hdaddy=gcf;
end

if(strcmp(get(hdaddy,'type'),'figure'))
    hfigkids=get(hdaddy,'children');
    haxes=[];
    for k=1:length(hfigkids)
        if(strcmp(get(hfigkids(k),'type'),'axes'))
            haxes=[haxes hfigkids(k)];
        end
    end
else
    haxes=hdaddy;
end

for kk=1:length(haxes)

    hkids=allchild(haxes(kk));

	for k=1:length(hkids)
		if(strcmp(get(hkids(k),'type'),'line')|strcmp(get(hkids(k),'type'),'hggroup'))
			lsize=get(hkids(k),'linewidth');
			set(hkids(k),'linewidth',linewidth*lsize);
            if(~strcmp(get(hkids(k),'type'),'hggroup'))
                msize=get(hkids(k),'markersize');
                set(hkids(k),'markersize',markersize*msize);
            end
            if(nargin==4)
                set(hkids(k),'color',kol)
            end
		end
	end
end
