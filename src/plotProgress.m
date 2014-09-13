function plotProgress(X,Y,Z,currentPlot,varargin)
shotIndex = [];
axisLimits = [];
cLimits = 'auto';

% check to see if figure exists
hf = findobj('Tag','SeimicExampleProgressPlot');
if isempty(hf)
    hf = figure('Tag','SeimicExampleProgressPlot');
    addpref('SeismicExample','figureHandle',hf);
    addpref('SeismicExample','timeHandle',[]);
    set(hf,'DeleteFcn',@(s,e)rmpref('SeismicExample'));
end
figure(hf)

% parse options
i = 1;
while i <= length(varargin)
    switch lower(varargin{i})
        case 'trace'
            shotIndex = varargin{i+1};
        case 'alim'
            axisLimits = varargin{i+1};
        case 'clim'
            cLimits = varargin{i+1};
        case 'time'
            % time to hours, minutes, seconds
            t = varargin{i+1};
            str = sec2hms(t);
            str = sprintf('Elapsed time: %s:%s:%s',str{1},str{2},str{3});
            ht = getpref('SeismicExample','timeHandle');
            if isempty(ht)
                % Create textbox for time
                ht = annotation(gcf,'textbox',[0.0 0.0 0.368 0.0385],...
                    'String',str,'FitBoxToText','on',...
                    'EdgeColor','none','Tag','ElapsedTimeLabel');
                setpref('SeismicExample','timeHandle',ht);
            else
                set(ht,'String',str);
            end
    end % case
    i = i+2;
end % for

% Plot
switch lower(currentPlot)
    case 'velocity'
        subplot(2,2,1)
        imagesc(X*1e-3,Y*1e-3,Z);
        %xlabel('Surface Distance (km)')
        ylabel('Depth (km)')
        title('Velocity')
        axis(axisLimits);
        caxis(cLimits);
        
    case 'stacked'
        subplot(2,2,2)
        imagesc(X*1e-3,Y*1e-3,Z);
        %xlabel('Surface Distance (km)')
        ylabel('Depth (km)')
        title('Migrated and Stacked Image')
        axis(axisLimits);
        caxis(cLimits);
        
    case 'migrated'
        subplot(2,2,4)
        imagesc(X*1e-3,Y,Z);
        %xlabel('Surface Distance (km)')
        ylabel('Time (s)')
        title('Migrated Shot')
        axis(axisLimits);
        caxis(cLimits);
        
    case 'shot'
        subplot(2,2,3)
        imagesc(X*1e-3,Y,Z);
        xlabel('Surface Distance (km)')
        ylabel('Time (s)')
        title(['Current Shot Record: ',num2str(shotIndex)])
        axis(axisLimits);
        caxis(cLimits);
        
        % add white window to velocity plot
        xrect = [X(1)  X(end) X(end) X(1) X(1)]*1e-3;
        hw = findobj('Tag','shotWindow');
        if isempty(hw)
            subplot(2,2,1), hold on
            z = axis;
            z = z(3:4);
            zrect = [z(end) z(end)  z(1)    z(1)  z(end)];
            plot(xrect,zrect,'w-','Tag','shotWindow')
            hold off
        else
            set(hw,'XData',xrect);
        end
end

function [str,hms] = sec2hms(t)
hms(1) = floor(t/3600);
hms(2) = floor((t-hms(1)*3600)/60);
hms(3) = floor(t-hms(1)*3600-hms(2)*60);

% covert to 00:00:00 format
str = cell(1,3);
for ii = 1:3
    if hms(ii) < 10
        str{ii} = ['0',num2str(hms(ii),'%2.0f')];
    else
        str{ii} = num2str(hms(ii),'%2.0f');
    end % if
end % for

