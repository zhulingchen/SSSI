function plotTrace(data)
% PLOTTRACE plot seismic data traces
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

[nSamples, nTraces] = size(data);

if (nSamples * nTraces == 0)
    fprintf('No traces to plot\n');
    return;
end

if (nSamples <= 1)
    fprintf('Only one sample per trace. Data set is not going to plot.\n');
    return;
end

times = (1:nSamples).';
location = 1:nTraces;

deflection = 1.25;

%% scale data
dx = (max(location)-min(location)) / (nTraces-1);
trace_max = max(abs(data));
scale = dx * deflection / (max(trace_max) + eps);
data = scale * data;

%% a figure window in portrait mode
hFig = figure;
set(hFig, 'Position', [288, 80, 900, 810], 'PaperPosition', [0.8, 0.5, 6.5, 8.0], ...
    'PaperOrientation', 'portrait', 'Color', 'w', 'InvertHardcopy', 'off');
figure(hFig)

set(gca, 'Position', [0.14, 0.085, 0.75, 0.79]);
set(hFig, 'Color', 'w');

axis([min(location)-deflection, max(location) + deflection, 1, nSamples]); hold on;
hCurAx = get(gcf, 'CurrentAxes');
set(hCurAx, 'ydir','reverse')
set(hCurAx, 'XAxisLocation','top');

%% plot data
ue_seismic_plot(times, data, location);
xlabel('Trace Number'); ylabel('Samples');
set(gca, 'XTick', location, 'XTickLabel', location);
grid on;
set(hCurAx, 'gridlinestyle', '-', 'box', 'on', 'xgrid', 'off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function plot seismic traces at horizontal locations controlled by location
function ue_seismic_plot(times, trace, location)

peak_fill = 'k';
trough_fill = '';

for ii = 1:size(trace, 2)
    
    y = trace(:, ii);
    
    chg = find(y(1:end-1).*y(2:end) < 0);
    x_zero = abs(y(chg) ./ (y(chg+1)-y(chg))) + times(chg);
    [x_data, idx] = sort([times(1); times; x_zero; times(end)]);
    y_data = [0; y; zeros(length(x_zero)+1, 1)];
    y_data = y_data(idx);
    
    h1 = fill(y_data(y_data >= 0) + location(ii), x_data(y_data >= 0), peak_fill);
    set(h1, 'EdgeColor', 'none');
    
    if ~isempty(trough_fill);
        h1 = fill(y_data(y_data <= 0) + location(ii), x_data(y_data <= 0), trough_fill);
        set(h1, 'EdgeColor', 'none');
    end
    
    
    plot([location(ii), location(ii)], [times(2), times(end)], 'w-')
    
    line(y_data(2:end-1) + location(ii), x_data(2:end-1), 'Color', 'k', ...
        'EraseMode', 'none', 'LineWidth', 0.5);
    
    
end