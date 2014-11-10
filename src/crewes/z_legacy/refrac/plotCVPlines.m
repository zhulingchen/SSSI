function plotCVPlines( action, axeslist, xi, xj, ti, tj )
[tmp numaxes] = size(axeslist);
k = 1;
if(isempty(xi))
   xi = NaN;
end
if(isempty(xj))
   xj = NaN;
end
linehandles = get(axeslist(1),'userdata');
[nlines tmp] = size(linehandles);
if( strcmp(action,'move') )
   thecolor = [0 0 0];
   if( ~isnan(xi) )
      thecolor = [0 0 1];    % Blue is the color of the 'i' lines
      newx = xi;
   elseif( ~isnan(xj) )
      thecolor = [0 1 0];    % Green is the color of the 'j' lines
      newx = xj;
   end
   
   str = sprintf('drawing lines at %f',newx);
   disp(str);
   for i=1:nlines
      l = linehandles(i);
      if( get(l,'color') == thecolor )
         set(l,'xdata',newx);
      end
   end
   drawnow;
end
if( strcmp(action,'draw') )
   drawmode = 'normal';
   if( nlines > 0 )
      delete(linehandles);
   end
   k = 1;
   for i=1:numaxes
      axes(axeslist(i));
      yr = get(gca, 'ylim');
      if( ~isnan(xi) )
         linehandles(k) = line([xi xi], yr, 'color', 'g' );
         set(linehandles(k),'erasemode','xor');
         k = k + 1;
      end
      if( ~isnan(xj) )
         linehandles(k) = line([xj xj], yr, 'color', 'b' );
         set(linehandles(k),'erasemode','xor');
         k = k + 1;
      end
   end
   axes(axeslist(1));
   if( ~isnan(xi) )
      linehandles(k)= line(xi, ti, 'linestyle', 'o', 'color', 'c');
      k = k+1;
   end
   if( ~isnan(xj) )
      linehandles(k) = line(xj, tj, 'linestyle', '+', 'color', 'c');
      k = k+1;
   end
   
   set(axeslist(1), 'userdata', linehandles);
end
