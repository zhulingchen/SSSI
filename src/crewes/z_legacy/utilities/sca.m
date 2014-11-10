function sca
% SCA: set current axis utility
%
% sca
%
% SCA is usefull with figures with multiple axes in setting the current axes as
% MATLAB seems to occaisionally screw this up

% search the figure for axes children
kids=get(gcf,'children');
axes=zeros(1,20); %allow for 20 in a figure
nax=0;
for k=1:length(kids)
   if(strcmp(get(kids(k),'type'),'axes'))
      nax=nax+1;
      axes(nax)=kids(k);
   end
end
p1=get(gca,'currentpoint');
% trap for matlab bug which does not set the current axes properly
if(nax>1)
   for k=1:nax
      ylim=get(axes(k),'ylim');
      xlim=get(axes(k),'xlim');
      p1=get(axes(k),'currentpoint');
      if( between(ylim(1),ylim(2),p1(1,2),2) & between(xlim(1),xlim(2),p1(1,1),2) )
	 set(gcf,'currentaxes',axes(k));
	 break;
      end
   end
end
return;

