function bhead = SEGY_ReadBinaryHeader(segy)

% bhead = SEGY_ReadBinaryHeader(segy)
% 'segy' is a SEGY structure, open for reading by SEG_OpenFile().
% bhead will be returned as a struct.
%
% Chad Hogan, 2008
%
% $Id: SEGY_ReadBinaryHeader.m,v 1.3 2009/07/24 15:49:19 kwhall Exp $

FILE = segy.FILE;

if(fseek(FILE, 3200, 'bof') ~= 0)
    % rewind the file to the binary header start
    disp(ferror(FILE));
    error('fseek failed in SEGY_ReadBinaryHeader() line 16');
end

bhead.jobid  = fread(FILE,1,   'int32');
bhead.lino  = fread(FILE,1,    'int32'); 
bhead.reno= fread(FILE,  1,    'int32');  
bhead.ntrpr = fread(FILE,  1,   'int16');
bhead.nart = fread(FILE,  1,    'int16');
bhead.hdt = fread(FILE,   1,    'uint16'); 
bhead.dto = fread(FILE,   1,    'uint16');
bhead.hns = fread(FILE, 1,      'uint16');
bhead.nso = fread(FILE, 1,      'uint16');
bhead.format = fread(FILE, 1,   'int16');
bhead.fold = fread(FILE, 1,     'int16');
bhead.tsort = fread(FILE, 1,    'int16');
bhead.vscode = fread(FILE, 1,   'int16');
bhead.hsfs = fread(FILE, 1,     'int16');
bhead.hsfe = fread(FILE, 1,     'int16');
bhead.hslen = fread(FILE, 1,    'int16');
bhead.hstyp = fread(FILE, 1,    'int16');
bhead.schn = fread(FILE, 1,     'int16');
bhead.hstas = fread(FILE,1,     'int16');
bhead.hstae = fread(FILE, 1,    'int16');
bhead.htatyp = fread(FILE, 1,   'int16');
bhead.hcorr = fread(FILE, 1,    'int16');
bhead.bgrcv = fread(FILE, 1,    'int16');
bhead.rcvm = fread(FILE, 1,     'int16');
bhead.mfeet = fread(FILE, 1,    'int16');
bhead.polyt = fread(FILE, 1,    'int16');
bhead.vpol = fread(FILE, 1,     'int16');

%skip unassigned hdr values
if(fseek(FILE, 3500, 'bof') ~= 0)
    % rewind the file to the binary header start
    disp(ferror(FILE));
    error('Error: fseek failed in SEGY_ReadBinaryHeader()');
end
bhead.rev = fread(FILE, 1,      'uint16');
bhead.trfix = fread(FILE, 1,    'uint16');
bhead.nthdr = fread(FILE, 1,    'int16');

%position file pointer at end of binary header
if(fseek(FILE, 3600, 'bof') ~= 0)
    % rewind the file to the binary header start
    disp(ferror(FILE));
    error('Error: fseek failed in SEGY_ReadBinaryHeader()');
end

%calculate number of bytes per trace
%assumes fixed trace length
bhead.tracebytelen = (1920 + (bhead.hns * 32)) / 8;

