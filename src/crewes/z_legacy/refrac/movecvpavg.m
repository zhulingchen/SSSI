% Function to allow manual editing of CVP average locations
%
% Usage   movecvpavg(action)
% where action is one of:
%
%  init  - Initialize function.  Turns on average editing.
%  fini  - Saves current CVP average locations in database (refdata)
% clear  - Clears button callback function assignments
function points = movecvpavg(action)
if(nargin<1 | strcmp(action, 'init'))    % set the button down function
  set(gcf,'windowbuttondownfcn','movecvpavg(''init'')');
end
% Get the data needed for this function
[shotrange shotINClist shot fullaxeslist stddev] = avgcvpinfo('get');
axeslist = fullaxeslist(1:5);
cvpavg = refdata('get', 'cvpavg');
shotcoord = refdata('get', 'shotcoord');
fbcoord = refdata('get', 'fbcoord');
fbtime = refdata('get', 'fbtime');
validx = find(~isnan(fbcoord(shot,:)));
xrange = fbcoord(shot,validx);
maxx = max(xrange);
minx = min(xrange);
shotx = shotcoord(shot);
% Get the location of the current mouse-click
pt=get(gca,'currentpoint');
userx = pt(1,1);
% Check to see if the user clicked in the 'i' shot side or the 'j' side.
iniside = 0;
if( userx > minx & userx < shotx )
  iniside = 1;
end
injside = 0;
if( userx < maxx & userx > shotx )
  injside = 1;
end
if( iniside == 0 & injside == 0 )
  disp('Selected point not in valid range.');
  return;
end
xi = cvpavg(shot,1);
xj = cvpavg(shot,2);
% Interpolate to find the time of each CVP point
ti = 0;
if( ~isnan(xi) )
  ti = interp1(fbcoord(shot,:),fbtime(shot,:),xi);
end
tj = 0;
if( ~isnan(xj) )
  tj = interp1(fbcoord(shot,:),fbtime(shot,:),xj);
end
  
if(strcmp(action,'init'))
  disp('Movecvpavg: initializing...');
  set(gcf,'windowbuttonmotionfcn','');
  set(gcf,'windowbuttonupfcn','');
  set(gcf,'windowbuttondownfcn','movecvpavg(''pick'')' );
end
if( strcmp(action,'pick') )
  % Selectiontype normal means button 1 - move CVP average to this point
  if(strcmp(get(gcf,'selectiontype'),'normal'))
    if( iniside )
      xi = userx;
      ti = interp1(fbcoord(shot,:),fbtime(shot,:),userx);
    else
      xj = userx;
      tj = interp1(fbcoord(shot,:),fbtime(shot,:),userx);
    end
    plotCVPlines( 'draw', axeslist, xi, xj, ti, tj ); 
  else
    str = sprintf('Placing CVP average at:',userx);
    disp(str);
  end
  
  set(gcf,'windowbuttonupfcn','movecvpavg(''fini'')');
  return;
end
if(strcmp(action,'fini'))
  str = sprintf('Movecvpavg: finishing shot %d ...', shot);  
  disp(str);
  if( iniside )
    xi = userx;
    ti = interp1( fbcoord(shot,:), fbtime(shot,:), userx );
    cvpavg(shot,1) = userx;
    refdata('set', 'cvpavg', cvpavg);
  elseif( injside )
    xj = userx;
    cvpavg(shot,2) = userx;
    refdata('set', 'cvpavg', cvpavg);
    tj = interp1(fbcoord(shot,:),fbtime(shot,:),userx);
  else
    disp('Point not in valid range.  Not set.');
  end
  plotCVPlines('draw',axeslist,xi,xj, ti, tj);
  
  set(gcf,'windowbuttonmotionfcn','');
  set(gcf,'windowbuttonupfcn','');
end
 
 
% if(strcmp(action,'save'))
%    objs = get(gca,'children');
%    n = 1;   % Counter for the number of picked points
%    [number x] = size(objs);
%    for i=1:number    % i is the counter for the number of axes children
%       if(strcmp(get(objs(i),'linestyle'),'*'))
%          point(n) = objs(i);
%          get(point(n),'xdata')   % This just prints X
%          n = n+1;
%       end
%    end
%    points = point;
% end
		
		
if( strcmp(action, 'clear'))
   set(gcf, 'windowbuttonmotionfcn', '');
   set(gcf, 'windowbuttonupfcn', '');
   set(gcf, 'windowbuttondownfcn', '');
end
