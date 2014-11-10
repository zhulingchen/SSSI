function SEGY_WriteStack(fname, stack, dt, separation, texttype, numtype)
% SEGY_WRITESTACK(fname, stack, dt, separation, texttype, numtype)
%
% Writes to the 'fname' file in SEGY, a 'stack' (NxM array, N samples 
% in M traces) with a sample interval of 'dt' seconds, a trace-to-trace 
% horizontal 'separation' distance, a 'texttype' equal to 'ebcdic' 
% or 'ascii', and a 'numtype' of 'l' for little-endian or 'b' for 
% big-endian.
%
% Chad Hogan, 2004
%
% $Id: SEGY_WriteStack.m,v 1.2 2004/07/30 21:23:35 kwhall Exp $
bhead  = SEGY_GetBinaryHeader;

if (strcmpi(texttype, 'ascii'))
    thead  = SEGY_GetTextHeader;
else
    thead = ascii2ebcdic(SEGY_GetTextHeader);
end

bhead.hdt = dt * 1000000; % dt is in s, but hdt in microseconds.
bhead.dto = bhead.hdt;

[tracesamps, numtraces] = size(stack);

% Now we set the number of samples.
bhead.hns = tracesamps;
bhead.nso = bhead.hns;

% Choose our endianness.
if (strcmpi(numtype, 'l'))
    machine = 'ieee-le';
    warning('Writing in non-standard little endian format');
else
    machine = 'ieee-be';
end

% open the file
[FILE, message] = fopen(fname, 'w', machine);
if (message); warning(message); end

SEGY_WriteTextHeader(FILE, thead);
segytrace = SEGY_GetTrace;
segytrace.id = 1;         % just one id in a file.
segytrace.dt = bhead.hdt; % time delta t
segytrace.ns = bhead.hns; % how many samples?

bhead.tsort = 4; % stacked
SEGY_WriteBinaryHeader(FILE, bhead);    
for i = 1:numtraces
    segytrace.cdp    = i; 
    segytrace.cdpt   = i; % trace number within gather
    segytrace.sx     = separation * (i - 1); % source location
    segytrace.gx     = segytrace.sx;         % group location
    segytrace.offset = 0;
    segytrace.data   = stack(:, i);
    SEGY_WriteTrace(FILE, segytrace, segytrace.ns);
end

fclose(FILE);
