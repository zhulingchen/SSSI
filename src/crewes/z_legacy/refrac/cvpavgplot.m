function cvpavgplot( graph, axishandle, xleft, yleft, xright, yright, xlimits )
axes(axishandle);
cla;
hold off;
if( length(yleft) > 0 )
   plot( xleft, yleft, 'y');
   hold on;
end
if( length(yright) > 0 )
   plot( xright, yright, 'c');
   hold on;
end
set( axishandle, 'xlim', xlimits );
hold off;
if( strcmp(graph, 'td') )            % time difference curves 
   string = 'Time difference (TD)';
   ylabel('Time (ms)');
end
if( strcmp(graph, 'median') )        % Plot the median filter
   set(axishandle, 'xtick', [] );
   string = 'Median filtered TD';
   ylabel('Time (ms)');
end
if( strcmp(graph, '1std') )          % Plot the 1st derivative
   set(axishandle, 'xtick', [] );
   string = '1st derivative of filtered TD';
end
if( strcmp(graph, '2ndd') )           % Plot the 2nd derivative
   set(axishandle,'xtick',[]);
   string = '2nd derivative of filtered TD';
end
text('units', 'normalized', 'position', [0.3 0.88], 'string', string);
