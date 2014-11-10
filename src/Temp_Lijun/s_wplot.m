function aux=s_wplot(seismic,varargin)
% Function plots seismic data in wiggle-trace form. 
% By default, positive deflections are filled black.
% Written by: E. R.: April 15, 2000
% Last updated: January 30, 2006: Add third option for cursor tracking
%
%             aux=s_wplot(seismic,varargin)
% INPUT
% seismic     seismic structure
% varargin    one or more cell arrays; the first element of each cell array 
%             is a keyword, the other elements are parameters. 
%             Presently, keywords are:
%     'annotation'  2-element or 3-element cell array. Defines a header 
%             mnemonic whose values are used to annotate the horizontal axis. 
%             Default: {'annotation','trace_no'}
%     'aindex annotation index; used to specify which traces to annotate with the
%             value of the header selected via keyword annotation.
%             EXAMPLE {'aindex',10:10:1000}   every tenth trace annotated, 
%                                             starting with trace 10
%             Default: {'aindex',[]}          Matlab-default annotation
%     'deflection'  trace deflection in units of trace spacing. 
%             Default: {'deflection',1.25}
%     'direction'  2-element cell array defining plot direction. Possible values are: 
%             left-to-right, 'l2r', and right-to-left, 'r2l'.  Default: {'direction','l2r')
%     'figure_only' specifies if only a figure window should be created; possible values:
%             'yes', 'no'. 
%             Default: {'figure_only','no'}
%     'figure'     Specifies if new figure should be created or if the seismic 
%             traces  should be plotted to an existing figure. Possible values
%             are 'new' and any other string. 
%             Default: {'figure','new'} 
%     'interpol'    2-element cell array which specifies how data should be interpolated in time.
%             Possible values are 'linear', 'cubic', and 'v5cubic'. 
%             'cubic'   - piecewise cubic Hermite interpolation
%             'v5cubic'  - the cubic interpolation from MATLAB 5, which does not
%                  extrapolate and uses 'spline' if X is not equally spaced.
%             Default: {'interpol','v5cubic'}
%     'orient' Plot orientation. Possible values are: 'portrait' and 'landscape'
%             Default: {'orient','landscape'} for more than 10 traces
%                      {'orient','portrait'}  for 10 or fewer traces
%     'peak_fill'   2-element cell array which specifies color of peak fill. Possible values are 
%             all permissible colors or the empty string. In the latter case peaks are not filled. 
%             Default: {'peak_fill','k'}
%     'pixels'      2-element cell array. Minimum number of interpolation points to use for display. 
%             The greater the number the smoother the display if 'interpol' is 'cubic'. 
%             Default: {'pixels',500}
%     'polarity'    2-element cell array. Possible values are 1 and -1;
%             Default: {'polarity',1}
%     'quality'     2-element cell array. Possible values are 'draft' and 'high'
%             'draft' quality (default) is faster and intended for screen and b/w hard copies; 
%             'high' quality is for color copies (there is no difference between
%             'draft' and 'high' quality for screen displays and b/w copies).
%             Default: {'quality','draft'}
%     'scale'       2-element cell array which specifies if traces should be scaled individually.
%             Possible values are 'yes', 'no', or the actual scale. This may be a scalar or a
%             vector whose number of elements is equal to the number of traces.
%             separation. The scale actually used must be obtained by specifying the output
%             argument "aux".
%             Default: {'scale','no'}
%     'spacing'     2-element cell array which specifies if traces should be   
%             equidistant ('equal') or non uniformly ('nonequal') spaced; 
%             in the latter case the header mnemonic used for annotation defines 
%             the trace-to-trace separation.
%             Default: {'spacing','equal'}
%     'times'       2-element or 3-element cell array 
%             {'times',vector of first and last time to plot} or ('times',first,last}. 
%             Default: {'times',seismic.first,seismic.last} which is
%                      equivalent to {'times',[seismic.first,seismic.last]}
%     'title'       2-element cell array consisting of the keyword 'title' and a title string;
%             no title is plotted if the title string is empty.
%             Default: {'title',seismic.name)
%     'traces'      2-element or 3-element cell array. The second element can be an array of 
%             trace numbers or it can be a string. If it is a string it can be a header 
%             mnemonic or it can contain a logical expression involving header values to 
%             include. A "pseudo-header" 'trace_no' can also be used.
%             If the second element is a string containing a header mnemonic there must 
%             be a third element containing a vector of values. (see "s_select")
%             Default:  {'traces',[]} which is equivalent to 
%                       {'traces',1:ntr} where ntr denotes the number of traces in the 
%                              input data set (ntr = size(seismic.traces,2))
%     'tracking' track cursor position; possible values are 'yes', 'no', and ''.
%             In the latter case a tracking button is created if the the
%             seismic is plotted in a new figure. Otherwise it is not.
%             Default: {'tracking',''}
%     'trough_fill' 2-element cell array which specifies color of trough fill. 
%             Possible values are all permissible colors or the empty string. 
%             In the latter case troughs are not filled. 
%             Default: {'trough_fill',''}
%     'wiggle'      2-element cell array which specifies color of wiggles. 
%             Possible values are all permissible colors or the empty string. 
%             In the latter case wiggles are not plotted. 
%             Default: {'wiggle','k'}
%
% OUTPUT
% aux        optional structure with scale factor(s) used (required if seismic
%            data are to be plotted with the same scaling in different plots)
%            'scale'   field with scale factor(s) used           
%            'figure_handle'  handle of the figure with the plot 
%
% EXAMPLE
%            seismic=s_data;
%            s_wplot(seismic,{'annotation','cdp'},{'direction','r2l'})

global S4M

if ~istype(seismic,'seismic')
   if isnumeric(seismic)
      seismic=s_convert(seismic,1,1);
      seismic.units='Samples';
   else
      error('First input argument must be seismic dataset or a matrix.')
   end
end

run_presets_if_needed

aux=[];

ntr=size(seismic.traces,2);

%     Set default values
param.annotation='trace_no';
if ntr < 11
   param.aindex=1:ntr;
elseif ntr < 21
   param.aindex=1:2:ntr;
else
   param.aindex=[];
end
param.deflection=1.25;
param.direction='l2r';
param.figure='new';
param.figure_only='no';
param.first=1;
param.inc=[];
param.interactive=1;
param.interpol='v5cubic';
param.orient=[];
param.peak_fill='k';
param.pixels=500;
param.polarity=1;
param.quality='draft';
param.scale='no';
param.spacing='equal';
param.times=[];
param.title=strrep(seismic.name,'_','\_');
param.traces=[];
param.tracking='';
param.trough_fill='';
param.wiggle='k';
param.wiggle_width=0.5;

%       Replace default values by input arguments
param=assign_input(param,varargin,'s_wplot');

history=S4M.history;  % Preserve value of global variable S4M.history

S4M.history=0;

if ~isempty(param.traces)
   if ~iscell(param.traces)
      seismic=s_select(seismic,{'traces',param.traces});
   else
      seismic=s_select(seismic,{'traces',param.traces{1},param.traces{2}});
   end
   ntr=size(seismic.traces,2);
end
if ~isempty(param.times)
   if length(param.times) <= 1
      seismic=s_select(seismic,{'times',param.times});
   elseif iscell(param.times)
      seismic=s_select(seismic,{'times',param.times{1},param.times{2}});
   else
      seismic=s_select(seismic,{'times',param.times(1),param.times(2)});
   end

end

if strcmp(param.peak_fill,'gray')
   param.peak_fill=[0.6,0.6,0.6];
elseif strcmp(param.peak_fill,'none')
   param.peak_fill=[];
end

if strcmp(param.trough_fill,'gray')
   param.trough_fill=[0.6,0.6,0.6];
elseif strcmp(param.trough_fill,'none')
   param.trough_fill=[];
end

if strcmp(param.wiggle,'gray')
   param.wiggle=[0.6,0.6,0.6];
elseif strcmp(param.wiggle,'none')
   param.wiggle=[];
end

S4M.history=history;

annotation=s_gh(seismic,param.annotation);
if ntr == 1
   uniform=logical(1);
   annotlog=logical(0);
else
   ddd=diff(annotation);
   annotlog=(max(ddd)==min(ddd)) & ddd(1) ~= 0 & isempty(param.aindex);
   if annotlog
      uniform=logical(0);
   else
      if strcmp(param.spacing,'equal')
         uniform=logical(1);        % Uniform trace-to-trace distance
      else
         uniform=logical(0);
      end
   end
end

[nsamp,ntr]=size(seismic.traces);
if nsamp*ntr == 0
   disp([' Alert from "s_wplot": Data set "',seismic.name,'" has no traces to plot'])
   return
end

%       Set indices of trace headers for annotation
if isempty(param.aindex)
   indices=1:max([round(ntr/7),1]):ntr;
else
   indices=param.aindex(param.aindex > 0 & param.aindex <= ntr);
end

if nsamp <= 1
   disp([' Only one sample per trace;  s_wplot did not plot data set "',seismic.name,'"'])
   return
end

%     Change polarity if requested
if param.polarity < 0
   seismic.traces=-seismic.traces;
end

%     Interpolate data if necessary
if (strcmpi(param.interpol,'v5cubic') | strcmpi(param.interpol,'cubic')) &  nsamp < param.pixels
   if isfield(seismic,'null') & isnan(seismic.null)
      for ii=1:ntr           % Replace last NaN at the top and first NaN at the bottom by zero
         idx=find(~isnan(seismic.traces(:,ii)));
         if ~isempty(idx)
            if idx(1) > 1
               seismic.traces(idx(1)-1,ii)=0;
            end
            if idx(end) < nsamp
               seismic.traces(idx(end)+1,ii)=0;
            end
         end
      end
   end
   npix=round(param.pixels/(nsamp-1))*(nsamp-1);
   dti=(seismic.last-seismic.first)/npix;
   times=(seismic.first:dti:seismic.last)';
   if isfield(seismic,'null')  &  S4M.matlab_version >= 7 % Turn off warnings caused by NaN's in seismic traces
      warning('off','MATLAB:interp1:NaNinY')
      yi=interp1(seismic.first:seismic.step:seismic.last,seismic.traces,times,param.interpol);
      warning('on','MATLAB:interp1:NaNinY')
   else
      yi=interp1(seismic.first:seismic.step:seismic.last,seismic.traces,times,param.interpol);
   end
   
else
%  dti=seismic.step;
   times=(seismic.first:seismic.step:seismic.last)';
   yi=seismic.traces;
end

trace_max=max(abs(yi));

%     Compute horizontal trace locations
if uniform
   xi=1:ntr;
else
   xi=annotation;
   if min(xi) == max(xi) & length(xi) > 1
      disp([' Header requested for annotation ("',param.annotation,'") is constant;'])
      disp('  please use uniform spacing')
      error(' Abnormal termination')
   end
end

%    Trace-to-trace separation
if ntr > 1
   dx=(max(xi)-min(xi))/(ntr-1);
else
   dx=1;
end

%     Scale data
if ischar(param.scale)
   if strcmpi(param.scale,'yes')
      scale=(dx*param.deflection)./(trace_max+eps)';
      yi=yi*spdiags(scale,0,ntr,ntr);
   else
      scale=dx*param.deflection/(max(trace_max)+eps);
      yi=yi*scale;
   end
else
  scale=param.scale;
  if length(scale) == 1
    yi=yi*scale;
  elseif length(scale) == ntr
    yi=yi*spdiags(scale(:),0,ntr,ntr);
  else
    error(' Scale factor must be scalar or vector of length equal to the number of traces')
  end
end

if nargout == 1
   aux.scale=scale;
end

%       Create figure window (unless an existing figure is to be used)
if strcmpi(param.figure,'new')
   if isempty(param.orient)
      if ntr > 10
         figure_handle=lfigure;
      else
         figure_handle=pfigure;
      end
   else
      if strcmpi(param.orient,'portrait')
         figure_handle=pfigure;
      else
         figure_handle=lfigure;
      end
   end
% 
   set(figure_handle,'Color','w')

   figure_export_menu(figure_handle);
   wseismic_scrollbar_menu(figure_handle,seismic,param.direction);

else
   figure_handle=gcf;
end

if nargout > 0
   aux.figure_handle=figure_handle;
end

set(gcf,'DoubleBuffer','on')

% dx_aver=(max(xi)-min(xi))/length(xi)
axis([min(xi)-dx*max([param.deflection,1]),max(xi)+dx*max([param.deflection,1]),seismic.first,seismic.last]);
hold on
ha=get(gcf,'CurrentAxes');
set(ha,'ydir','reverse')
set(ha,'XAxisLocation','top');

%    Handle reversal of plot direction
if strcmpi(param.direction,'r2l')
   set(ha,'xdir','reverse')
   yi=-yi;
   temp=param.trough_fill;
   param.trough_fill=param.peak_fill;
   param.peak_fill=temp;
elseif ~strcmpi(param.direction,'l2r')
   error(['Keyword for plot direction is wrong (',param.direction,')'])
end

if strcmpi(param.figure_only,'yes')
   return
end

%    Check if there are NaNs in the seismic data
if isfield(seismic,'null')
   isnull=1;
else
   isnull=0;
end


%    Plot data
if strcmpi(param.quality,'draft')
   ue_seismic_plot(times,yi,xi,isnull,param) 
else
   ue_seismic_plot_ps(times,yi,xi,isnull,param)
end

if (isempty(param.tracking) & strcmpi(param.figure,'new'))  |  ...
    strcmpi(param.tracking,'yes')

%	Add cursor tracking
   [dummy,xinfo]=s_gh(seismic,param.annotation);   %#ok First output argument is not required
   y=linspace(seismic.first,seismic.last,nsamp);
   yinfo=info4time(seismic);
   if ~uniform
      xi2use=xi;
   else
      xi2use=s_gh(seismic,param.annotation);
   end
   initiate_3d_tracking4seismic(seismic.traces,xi2use,y,xinfo,yinfo,{'amplitude','','Amplitude'})
else
   yinfo=info4time(seismic);
   ylabel([yinfo{3},' (',yinfo{2},')'])
end

%	Add annotation of horizontal axis
if ~annotlog
   v=axis;
   xil=xi(indices);
   annol=annotation(indices);

   add_xaxis([xil',annol'],'',{'location','top'});
end

%	Add label of horizontal axis
%       (if 5 or more traces or annotation other than trace number)
if ~strcmpi(param.annotation,'trace_no') | ntr > 4
    xtext=s_gd(seismic,param.annotation);
    hunits=s_gu(seismic,param.annotation);

    if ~isempty(xtext) & ~strcmpi(xtext,'not available')
       if ~isempty(hunits) & ~strcmpi(hunits,'n/a')
          xtext=[xtext,' (',hunits,')'];
       end
       hxlabel=xlabel(xtext);
       v=get(hxlabel,'Position'); % Preserve position of x-axis label
    else
       hxlabel=xlabel('');
    end
else
    hxlabel=xlabel('');
    set(gca,'XTickLabel',[])
end


if uniform & (~strcmpi(param.annotation,'trace_no') | ntr > 4)
%  set(gca,'XTickLabel',[]);	% Delete original annotation but
   try                          % keep x-axis label at original position
      set(hxlabel,'Position',v);    
   catch
   end
end


grid on; 
set(ha,'gridlinestyle','-','box','on','xgrid','off')

%  Title
if ~isempty(param.title)
   if iscell(param.title)		% Handle multi-line titles
      mytitle(param.title{1})
   else
      mytitle(param.title)
   end
end

if strcmpi(param.figure,'new')
   timeStamp
end

zoom off
hold off
   
if nargout == 0
   clear aux
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ue_seismic_plot(times,traces,location,isnull,param)
% Function plot seismic traces at horizontal locations controlled by location

% global param
if ~isempty(param.wiggle)
   wiggle=1;
else
   wiggle=0;
end

dt=times(2)-times(1);
times_orig=times;

for ii=1:size(traces,2)
   skip=0;        
   if isnull
      idx=find(~isnan(traces(:,ii)));
      if isempty(idx)
         skip=1;
      else
         y=traces(idx,ii);
         times=times_orig(idx);
      end
   else
     y=traces(:,ii);
   end

              if ~skip
   chg=find(y(1:end-1).*y(2:end) < 0);
   x_zero=abs(y(chg)./(y(chg+1)-y(chg)))*dt+times(chg);
   [x_data,idx]=sort([times(1);times;x_zero;times(end)]);
   y_data=[0;y;zeros(length(x_zero)+1,1)];
   y_data=y_data(idx);

   if ~isempty(param.peak_fill);
      pos=find(y_data >= 0);
      h1=fill(y_data(pos)+location(ii),x_data(pos),param.peak_fill); 
      set(h1,'EdgeColor','none')
   end

   if ~isempty(param.trough_fill);
      neg=find(y_data <= 0);
      h1=fill(y_data(neg)+location(ii),x_data(neg),param.trough_fill); 
      set(h1,'EdgeColor','none')
   end

   if ~isempty(param.trough_fill) | ~isempty(param.peak_fill)
      plot([location(ii),location(ii)],[times(2),times(end)],'w-')
   end
   if wiggle
      line(y_data(2:end-1)+location(ii),x_data(2:end-1),'Color',param.wiggle, ...
           'EraseMode','none','LineWidth',param.wiggle_width);
   end
              end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ue_seismic_plot_ps(times,traces,location,isnull,param)
% Function plot seismic traces at horizontal locations controlled by location
% Avoids vertical lines from top to bottom on some PS interpreters
% Much slower than ue_seismic_plot

% global param
if ~isempty(param.wiggle)
  wiggle=1; 
else
  wiggle=0;
end 

dt=times(2)-times(1);
times_orig=times;

for ii=1:size(traces,2)
  skip=0;        
  if isnull
    idx=find(~isnan(traces(:,ii)));
    if isempty(idx)
      skip=1;
    else
      y=traces(idx,ii);
      times=times_orig(idx);
    end
  else
    y=traces(:,ii);
  end
              if ~skip
  chg=find(y(1:end-1).*y(2:end) < 0);
  x_zero=abs(y(chg)./(y(chg+1)-y(chg))).*dt+times(chg);
  [x_data,idx]=sort([times(1);times;x_zero;times(end)]);
  y_data=[0;y;zeros(length(x_zero),1);0];
  y_data=y_data(idx);

  if ~isempty(param.peak_fill);
    pos=find(y_data >= 0);  
%    y_pos=y_data(pos); pos=find(y_data >= 0);
    ind1=find(diff(pos) > 1); 
    ind=[pos(1),reshape([pos(ind1)';pos(ind1+1)'],1,2*length(ind1)),pos(end)];
%    x_pos=x_data(pos);
    for kk=1:2:length(ind)-1;
      index=ind(kk):ind(kk+1);
      fill(y_data(index)+location(ii),x_data(index),param.peak_fill, ...
             'erasemode','none','EdgeColor','none'); 
    end
  end
  if ~isempty(param.trough_fill);
    neg=find(y_data <= 0);
    ind=find(diff(neg) > 1); 
    ind=[neg(1),reshape([neg(ind)';neg(ind+1)'],1,2*length(ind)),neg(end)];
    for kk=1:2:length(ind)-1
      index=ind(kk):ind(kk+1);
      fill(y_data(index)+location(ii),x_data(index),param.trough_fill, ...
             'erasemode','none','EdgeColor','none'); 
    end
  end

  if wiggle
     line(y_data(2:end-1)+location(ii),x_data(2:end-1),'Color',param.wiggle, ...
         'erasemode','none','LineWidth',param.wiggle_width); 
  end
              end
end
box on
