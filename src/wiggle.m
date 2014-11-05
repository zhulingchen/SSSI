%WIGGLE Display data as wiggles.
%   WIGGLE(C) displays matrix C as wiggles plus filled lobes, which is a
%   common display for seismic data or any oscillatory data. A WIGGLE
%   display is similar to WATERFALL, except that the Z heights are
%   projected onto the horizontal plane, meaning that a WIGGLE display is
%   actually a 2D display.
%
%   WIGGLE(C) plots each column of C as curves vertically. How much the
%   deviation is from the jth vertical is given by ith row. It C(i,j) is
%   positive then the curve bends to the right, otherwise it bends to the
%   left. Best visual effect is achieved when each column has no DC
%   component, i.e., when its sum is zero, or at least close to it.
%
%   WIGGLE(X,Y,C), where X and Y are vectors, rescale the axes to match X
%   and Y values. When WIGGLE has one or three input, its usage is very 
%   similar to IMAGE or IMAGESC, except for WIGGLE's properties, described 
%   below.
%
%   WIGGLE(...,S) allows some control of wiggles properties, in which S is
%   a string with up to six characters, each one setting up a property,
%   which can be Extrusion, Polarity, Direction, Wiggle Color and Lobes 
%   Color.
%
%   Extrusion: '0', '1', ..., '9'
%     It is how much a wiggle overlaps its neighbor. The default, '0', 
%     means there is no overlap. For instance, E='1' means that a
%     wiggle overlaps only one neighbor and so on. Observe that only 
%     positve numbers are allowed.
%
%   Polarity: '+' or '-'
%     It defines which side of the wiggle is going to be filled.
%     Polarity '+' (default) means that the positive lobes are filled,
%     while '-' means that the negative lobes are filled.
%
%   Wiggles Direction: 'v' or ' h'
%     It specifies if the wiggles are vertical ('v'), which is the default,
%     or horizontal ('h').
%
%   Wiggles Color: 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w' or 'i'
%     It defines the wiggle color. This property uses the same color key
%     as in function PLOT. The default is black, i.e., 'k'. In order to
%     suppress the wiggle, the character 'o' must be used. Options 'i',
%     invisible, means that no wiggle is displayed.
%
%   Lobes color: 'B', 'G', 'R', 'C', 'M', 'Y', 'K', 'W', 'X', '*' or 'I'
%     It defines filling color for the right and left lobes. This property 
%     uses mostly the same color key as used in command PLOT, but in 
%     uppercase. 
%     The order of preference is right before left lobes, meaning that if 
%     only one character is given, the right lobes are painted with the
%     corresponding color and the left lobes are invisible. In this way,
%     the default is left lobes invisible and right lobes in black, i.e., 
%     the string contains 'K' or 'IK'. Notice that a single 'I' implies 
%     that both lobes are invisible, because by default the left lobes are
%     already invisible. If two capital characters are given, the first 
%     one defines the color of the left lobes and the second one of the 
%     right lobes. 
%     There are still two special coloring options which use the current 
%     figure colormap: '*' means that the lobes ares filled with a variable 
%     color that is horizontally layered and 'X' means that the lobes are 
%     filled with a variable color that is vertically layered.
%
%   Note that the characters that build string S can be given in any
%   order, but respecting the range and precedence of each propetry:
%
%   Extrusion       Polarity         Direction               Color    
%                                                     Wiggles       Lobes 
%   0 - (default)   + - right lobe   v - vertical     b     blue        B
%   1               - - left lobe    h - horizontal   g     green       G
%   2                                                 r     red         R
%   ...                                               c     cyan        C
%                                                     m     magenta     M
%                                                     y     yellow      Y
%                                                     k     black       K
%                                                     w     white       W
%                                                     i     invisible   I
%                                                           colormap_h  X
%                                                           colormap_v  *
%
%   Examples
%      C = 2*rand(200,20) - 1;
%      D = filter([1,1,1],1,C);
%      figure(1); clf
%      subplot(131); wiggle(D);          title('default: right lobes in black')
%      subplot(132); wiggle(D, '-');     title('left lobes in black')
%      subplot(133); wiggle(D, 'yBR');   title('magenta lines, blue and red lobes')
%      figure(2); clf
%      subplot(311); wiggle(D', 'h');    title('Horizontal wiggles')
%      subplot(312); wiggle(D', 'hB');   title('Horizontal wiggles with blue lobes')
%      subplot(313); wiggle(D', 'hbMR'); title('Hor. wiggles, blue lines with magenta and red lobes')     
%      figure(3); clf
%      subplot(131); wiggle(D, 'BR');    title('Blue and red lobes')                  
%      subplot(132); wiggle(D, 'XX');    title('Horizontally filled colormapped lobes')                
%      subplot(133); wiggle(D, '1.5**'); title('Vertically filled colormapped lobes with 1.5 of extrusion') 
%
%   See also IMAGESC, IMAGE, WATERFALL, PLOT, LINESPEC

% by Rodrigo S. Portugal (rosoport.matlab@gmail.com)
% last revision: 18/10/2012

% TODO:
%   1) Implement the invisible option (DONE)
%   2) Implement the variable color option (DONE)
%   3) Implement wiggles overlaid imagesc (DONE)
%   3) Pair-wise options for fine tuning of properties (NOT YET)
%   3) Inspect the code (DONE)
%   4) Optimize (+/-)
%   5) Test examples (DONE)


function wiggle(varargin)

switch (nargin)
    
    case 0
        error('Too few input arguments.');
        
        % Only data
    case 1
        data = check_data(varargin{1});
        [nr,nc] = size(data);
        xr = 1:nr;
        xc = 1:nc;
        prop = '';
        
        % Data and properties
    case 2
        data = check_data(varargin{1});
        [nr,nc] = size(data);
        xr = 1:nr;
        xc = 1:nc;
        prop = check_properties(varargin{2});
        
        % Domain vectors and data
    case 3
        xr   = check_data(varargin{1});
        xc   = check_data(varargin{2});
        data = check_data(varargin{3});
        prop = '';
        
        % Domain vectors, data and properties
    case 4
        xr   = check_data(varargin{1});
        xc   = check_data(varargin{2});
        data = check_data(varargin{3});
        prop = check_properties(varargin{4});
        
    otherwise
        error('Too many input arguments.');
        
end

[extr, pol, dir, wcol, vcoll, vcolr] = extract_properties(prop);

wiggleplot(data, xr, xc, extr, wcol, vcoll, vcolr, dir, pol)

end

%--------------------------------------------------------------------------
function data = check_data(v)

if isnumeric(v),
    data = v;
else
    error('Wiggle: value must be numeric');
end

end

%--------------------------------------------------------------------------
function prop = check_properties(v)

if ischar(v)
    prop = v;
else
    error('Wiggle: properties must be a string.');
end

end

%--------------------------------------------------------------------------
function check_rgb(color)

   if any(isnan(color))  || any(color > 1) || ...
      size(color,1) ~= 1 || size(color,2) ~= 3 
         error(['Wiggle: color must be a numeric 1 x 3 matrix,  '...
                'following RGB system'])
   end
   
end

%--------------------------------------------------------------------------
function [extr, pol, dir, wcol, vcoll, vcolr] = extract_properties(prop)


wcol = 'k';
vcol = 'k';
dir  = 0;
pol  = 1;

extr  = extractfloat(prop);
if isnan(extr), extr = 1.0; end

indv = 1;

for ip = 1:length(prop)
    
    p = prop(ip);
    
   
    if p == '+' || p == '-'
        if p == '+' 
            pol = 1;
        else
            pol = -1;
        end
        continue
    end
    
    if p == 'h' || p == 'v' || p == 'H' || p == 'V'
        if p == 'v' || p == 'V',
            dir = 0;
        else
            dir = 1;
        end
        continue
    end
    
    if p == 'b' || p == 'g' || p == 'r' || p == 'c' || p == 'm' || ...
       p == 'y' || p == 'k' || p == 'w' || p == 'i' 
        wcol = p;
        continue
    end
    
    if p == 'B' || p == 'G' || p == 'R' || p == 'C' || p == 'M' || ...
       p == 'Y' || p == 'K' || p == 'W' || p == 'I' || p == 'X' || p == '*' 
        vcol(indv) = lower(p);
        indv = 2;
        continue
    end
    
end

wcol  = str2rgb(wcol);

if length(vcol) == 2
    vcoll = str2rgb(vcol(1));
    vcolr = str2rgb(vcol(2));
else
    vcoll = NaN;
    vcolr = str2rgb(vcol(1));
end


end

%--------------------------------------------------------------------------

function f = extractfloat(str)

fs = blanks(length(str));
foundpoint = false;
count = 1;
for i = 1: length(str)
    
    cs = str(i);
    if cs == '.' && ~foundpoint
        fs(count) = cs;
        count = count + 1;
        foundpoint = true;
        continue
    end
    
    c  = str2double(cs);
    if ~isnan(c) && isreal(c)
        fs(count) = cs;
        count = count + 1;
    end

end

f = str2double(fs);

end

%--------------------------------------------------------------------------

function rgb = str2rgb(s)

switch s
    
    case 'r'
        rgb = [1, 0, 0];
    case 'g'
        rgb = [0, 1, 0];
    case 'b'
        rgb = [0, 0, 1];
    case 'c'
        rgb = [0, 1, 1];
    case 'm'
        rgb = [1, 0, 1];
    case 'y'
        rgb = [1, 1, 0];
    case 'k'
        rgb = [0, 0, 0];
    case 'w'
        rgb = [1, 1, 1];
    case 'i'
        rgb = NaN;
    case 'x'
        rgb = 1;
    case '*'
        rgb = 2;
    otherwise
        rgb = NaN;
end
end

%--------------------------------------------------------------------------

% WIGGLEPLOT plots seismic data as wiggles and variable area
%  WIGGLEPLOT(data, y, x, k, wc, vcl, vcr, d, p)
%  INPUT  DESCRIPTION                              TYPE      SIZE        
%   data  oscilatory data                          matrix    (nx x nt)
%   y     vertical coordinates                     vector    ( 1 x ny)   
%   x     horizontal coordinates                   vector    ( 1 x nx)   
%   k     extrusion                                scalar    ( 1 x 1 )   
%   wc    wiggle color                             matrix    ( 1 x 3 )   
%         wc=0 supress the wiggle
%   vcl   variable area color (left lobe)          matrix    ( 1 x 3 )   
%         vcl=0 or NaN suppress the left variable area
%   vcr   variable area color                      matrix    ( 1 x 3 )   
%         vcr=0 or NaN suppress the right variable area
%   d     wiggle direction                         scalar    ( 1 x 1 )   
%         d=0 for vertical
%         d=1 for horizontal
%   p     polarity                                 scalar    ( 1 x 1 )   
%         p=1  (SEG convention) positive area
%         p=-1 (EAGE convention) negative area
%
%   See also IMAGESC, IMAGE

function wiggleplot(data, xr, xc, extr, wcol, vcoll, vcolr, dir, pol)


[nrows,ncols] = size(data);


if ncols ~= length(xc) || nrows ~= length(xr)
    error('Input data must have compatible dimensions');
end
if dir == 1
    data = data';
    extr = -extr;
    aux  = xc;
    xc    = xr;
    xr    = aux;
end
nxr = length(xr);
nxc = length(xc);

if nxc > 1,
    dx = xc(2)-xc(1);
else
    dx = 1;
end

% ir = nxr-1:-1:2;

plot_var = true;
plot_val = true;

if length(vcoll) == 1
    if isnan(vcoll) || vcoll == 0
        plot_val = false;
    elseif vcoll ~= 1 && vcoll ~= 2 
        error('wiggleplot: color not recognized.')
    end
else
    check_rgb(vcoll)
end

if length(vcolr) == 1
    if isnan(vcolr) || vcolr == 0
        plot_var = false;
    elseif vcolr ~= 1 && vcolr ~= 2 
        error('wiggleplot: color not recognized.')
    end
else
    check_rgb(vcolr)
end

plot_w  = true;

if length(wcol) == 1
    if isnan(wcol) || wcol == 0
        plot_w = false;
    else
        error('wiggleplot: color not recognized.')
    end
else
    check_rgb(wcol)
end

if pol == -1,
    data = -data;
end

count = make_counter(nxc, dir);

[plotdir, patchdirl, patchdirr] = select_functions(dir, plot_w, plot_val, plot_var);

scale = pol * extr * dx / (max(data(:))+eps);

hold on

xwig = zeros(nxc,2*nxr+1); 
ywig = zeros(nxc,2*nxr+1);
if vcoll == 2
    chunk = zeros(nxc,4*nxr+3);
else
    chunk = zeros(nxc,2*nxr+2);
end
xval = chunk; yval = chunk; %cval = chunk;
xvar = chunk; yvar = chunk; %cvar = chunk;

for ic = count(1):count(2):count(3)
    
    
    [xr1, d, sl1, sr1, ~, ~] = make_lines(xr, data(:,ic)');
    xr2 = xr1;
    sr2 = sr1;
    sl2 = sl1;
    
    if vcolr == 1
        colr = [0, d, 0];
    elseif vcolr == 2
        colr = [0, d, 0, fliplr(d), 0];
        xr2  = [xr1, fliplr(xr1(1:end-1))];
        sr2   = [sr1, zeros(1,length(sr1)-1)];
    else
        colr = vcolr;
    end
    
    xvar(ic,:) = xc(ic) + scale * sr2;
    yvar(ic,:) = xr2;
    cvar(ic,:) = colr;
    
    
    if vcoll == 1
        coll = [0, d, 0];
    elseif vcoll == 2
        coll = [0, d, 0, fliplr(d), 0];
        xr2  = [xr1, fliplr(xr1(1:end-1))];
        sl2   = [sl1,  zeros(1,length(sl1)-1)];
    else
        coll = vcoll;
    end
    
    xval(ic,:) = xc(ic) + scale * sl2;
    yval(ic,:) = xr2;
    cval(ic,:) = coll;
    
    
    xwig(ic,:) = [xc(ic) + scale * d, NaN];
    ywig(ic,:) = [xr1(2:end-1), NaN];
    
    
end

patchdirl(xval', yval', cval');
patchdirr(xvar', yvar', cvar');
plotdir(xwig,ywig,wcol);

hold off

if dir == 1,
    axis([min(xr), max(xr), min(xc) - extr * dx, max(xc) + extr * dx])
else
    axis([min(xc) - extr * dx, max(xc) + extr * dx, min(xr), max(xr)])
end

set(gca,'NextPlot','replace', 'XDir','normal','YDir','reverse')

end

%--------------------------------------------------------------------------
function counter = make_counter(nxc, dir)

if dir == 1,
    ic_first = 1;
    ic_last  = nxc;
    ic_step  = 1;
else
    
    ic_first = nxc;
    ic_last  = 1;
    ic_step  = -1;
end

counter = [ic_first, ic_step, ic_last];

end

%--------------------------------------------------------------------------
function [fw, fval, fvar] = select_functions(dir, plot_w, plot_val, plot_var)

fw  = @plot_vert;
fval = @patch_vert;
fvar = @patch_vert;

if dir == 1,
    fw = @plot_hor;
    fval = @patch_hor;
    fvar = @patch_hor;
end

if plot_w == false
    fw = @donothing;
end

if plot_val == false
    fval = @donothing;
end

if plot_var == false
    fvar = @donothing;
end
    
if (plot_val == false) && (plot_var == false) && (plot_w == false),
    disp('wiggle warning: neither variable area and wiggles are displayed');
end

end


%--------------------------------------------------------------------------

function [t1, d, sl, sr, ul, ur] = make_lines(t, s)

nt = length(t);
dt = t(2) - t(1);

ds = s(2:nt)-s(1:nt-1);

r = ((s(2:nt)<0) .* (s(1:nt-1)>0)) + ((s(2:nt)>0) .* (s(1:nt-1)<0));
a = r .* s(2:nt)./(ds+eps) + (1-r).*0.5;

tc(1:2:2*nt-1) = t;
tc(2:2:2*nt-2) = t(2:nt) - dt * a;
sc(1:2:2*nt-1) = s;
sc(2:2:2*nt-2) = 0.5*(s(1:nt-1)+ s(2:nt)).*(1-r);

t1 = [t(1),    tc, t(nt), t(nt)];
d0 = [0,       sc, s(nt), 0];
dr = [max(sc), sc, s(nt), max(sc)];
dl = [min(sc), sc, s(nt), min(sc)];


sl = min(0,d0);
sr = max(0,d0);

ul = min(0,dl);
ur = max(0,dr);

d = d0(2:end-1);
end

%--------------------------------------------------------------------------

function donothing(~,~,~)
end


%--------------------------------------------------------------------------

function ph = plot_vert(x,y,col)
x = reshape(x',1,numel(x));
y = reshape(y',1,numel(y));
ph = line(x,y,'linewidth',0.25,'color', col);
end

%--------------------------------------------------------------------------

function ph = plot_hor(x,y,col)
x = reshape(x',1,numel(x));
y = reshape(y',1,numel(y));
ph = line(y,x,'linewidth',0.25,'color', col);
end


%--------------------------------------------------------------------------

function fh = patch_vert(x,y,col)
if size(col,1) == 3,
    col = col(:,1)';
end
fh = patch(x,y,col,'edgecolor','none');
hold on
end

%--------------------------------------------------------------------------

function fh = patch_hor(x,y,col)
if size(col,1) == 3,
    col = col(:,1)';
end
fh = patch(y,x,col,'edgecolor','none');
hold on
end



