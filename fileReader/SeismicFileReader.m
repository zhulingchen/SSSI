classdef SeismicFileReader < handle
%SeismicFileReader Seismic file reader class.
%
%   SeismicFileReader is a handle class that provides functionality for
%   reading a binary seismic file.  Seismic files generally contain file
%   headers that describe the file contents and seiscmic traces.  The
%   seismic traces typically contain a trace header followed by the trace
%   data.  This is a generic class that can be used to construct file
%   readers for different formats (e.g. SEG2, SEGY, SEGD).
%
%   SeismicFileReader properties:
%       FileName            - name of the seismic binary file
%       FilePath            - path of the file
%       FileFormat          - file format description
%       FileHeader          - Header information from the seismic file
%       MachineFormat       - Binary file format, big/little endian
%       NumberOfTraces      - Number of traces in the file
%       TraceDescription    - Trace information from the seismic file
%
%  SeismicFileReader methods:
%      SeismicFileReader    - constructor, create a SeismicFileReader
%                             object
%      readFileHeader       - read in mixed format file header
%      readTraceHeader      - read in mixed format trace header
%      readTraceData        - read in numeric trace data
%
%   See also SEG2FILEREADER, SEGYFILEREADER
   
    properties (SetAccess = protected) 
%FILENAME Name of the binary file.
%   The FILENAME property records the file name and is used for opening the
%   file for reading.

        FileName = '';
        
%FILEPATH Path of the binary file.
%   The FILEPATH property contians the file location as a string.

        FilePath = '';
%FILESIZE size in bytes of file.
%   The FILESIZE property contians the file size in bytes.

        FileSize = '';
     
%FILEFORMAT File format information.
%   Te FILEFORMAT property is used to provide a descriptions of the file
%   format (e.g. seg2 version 0)
        
        FileFormat = '';

%FILEHEADER File descriptive information.
%   The FILEHEADER property contains the information contained in the file
%   header.

        FileHeader = {};
        
%MACHINEFORMAT Format of the binary file.
%   The MACHINEFORMAT property contains information about the binary file
%   type, big-endian or little-endian.

        MachineFormat = 'native';
        
%NUMBEROFTRACES Number of traces in the file.
%   The NUMBEROFTRACES property reports number of traces in the file, an
%   integer value.

        NumberOfTraces = NaN;

%TRACEDESCRIPTION Trace information.
%   The TRACEDESCRIPTION property contains information about the trace data
%   inf the file.

        TraceDescription = {};
        
    end % public properties
    
    % Property data is private to the class
    properties (SetAccess = private, GetAccess = private)
        FileID;                   % file id for reading
    end % private properties
          
    % Public Class Methods
    methods
        function obj = SeismicFileReader(varargin)
%SEISMICFILEREADER creates a file reader object.
%   SEISMICFILEREADER creates an empty SEISMICFILEREADER object.
%
%   SEISMICFILEREADER(FILENAME) creates a SEISMICFILEREADER by opening and
%   FILENAME for reading
%
%   See also READFILEHEADER, READTRACEHEADER, READTRACEDATA
            switch nargin
                case 0 % create object but don't link to a file
                    return; % return default object
                case 1 % filename and path provided
                    fileName = varargin{1};
                otherwise
                    error('SeismicFileReader:TooManyInputs',...
                        'Input should be a filename or empty.')
            end
            openFile(obj,fileName)
        end % constructor
        
        function h = readFileHeader(obj,fmt,varargin) 
%READFILEHEADER reads the header fields.
%   H = READFILEHEADER(OBJ,FMT,OFFSET) reads in the file header
%   information for OBJ (the SEISMICFILEREADER object) according to the
%   format specified in FMT.  It returns the information in a structure H
%   that contains filednames defined by FMT.  FMT is a cell array that
%   contains the binary format (e.g. uint8) and the variable names of the
%   data {FORMAT,{VARNAME1, VARNAME2,...},...} or for creating an vector
%   using a single variable name {FORMAT,{NUMBER,VARNAME},...}.  OFFSET is
%   the byte offset from the beginning of file to begin the read from.
%
%   Example:
%   % read in the header for a SEG2 formatted file (partial read)
%   s = SeismicFileReader('WFLT0001')
%   fmt = {'uint16',{'FileDescriptor','RevisionNumber',...
%                     'TracePointerSize','NumberOfTraces'};
%          'uint8', {'StringTerminatorSize',...
%                    'FirstStringTerminator',...
%                    'SecondStringTerminator',...
%                    'SizeofLineTerminator',...
%                    'FirstLineTerminator',...
%                    'SecondLineTerminator' };
%           };
%   fileHeader = readFileHeader(s,fmt)
%
%   % read in the file header information and lump terminator information
%   % into a single array named Terminators
%   fmt = {'uint16',{'FileDescriptor','RevisionNumber',...
%                     'TracePointerSize','NumberOfTraces'};
%          'uint8', {6, 'Terminators'} };
%   fileHeader = readFileHeader(s,fmt)
%
%   See also SEISMICFILEREADER, READTRACEHEADER, READTRACEDATA
            
            h = readHeader(obj,fmt,varargin{:});
            
        end % readFileHeader
        
        function h = readTraceHeader(obj,fmt,varargin)
%READTRACEHEADER reads the header fields.
%   H = READTRACEHEADER(OBJ,FMT,OFFSET) reads in the file header
%   information for OBJ (the SEISMICFILEREADER object) according to the
%   format specified in FMT.  It returns the information in a structure H
%   that contains filednames defined by FMT.  FMT is a cell array that
%   contains the binary format (e.g. uint8) and the variable names of the
%   data {FORMAT,{VARNAME1, VARNAME2,...},...} or for creating an vector
%   using a single variable name {FORMAT,{NUMBER,VARNAME},...}.  OFFSET is
%   the byte offset from the beginning of file to begin the read from.
%
%   Example:
%   % read in part of the trace header for a SEG2 formatted file.
%   s = SeismicFileReader('WFLT0001')
%
%   % read in the file header information to get the pointers to the trace
%   % data
%   fmt = {'uint16',{'FileDescriptor','RevisionNumber',...
%                     'TracePointerSize','NumberOfTraces'};
%          'uint8', {6, 'Terminators'} };
%   fileHeader = readFileHeader(s,fmt)
%
%   % get the number of traces and the trace locations in the file
%   fmt = {'uint32', {fileHeader.NumberOfTraces, 'TracePointers'}};
%   offset = 32; % bytes from beginning of file
%   locInFile = readFileHeader(s,fmt,offset);
%
%   % read in the header for trace 1
%   fmt = {'uint16',{'TraceDescriptor','TraceBlockSize'};...
%          'uint32',{'TraceDataSize','NumberOfSamples'};...
%          'uint8', {'TraceDataFormatCode'} };
%   offset = locInFile.TracePointers(1);      
%   traceHeader = readTraceHeader(s,fmt,offset)
%
%   See also SEISMICFILEREADER, READFILEHEADER, READTRACEDATA
            h = readHeader(obj,fmt,varargin{:});
        end % read Trace Header
        
        function d = readTraceData(obj,fmt,varargin)
%READTRACEDATA reads the trace numeric data
%   H = READTRACEDATA(OBJ,TRACES,OFFSET) reads in the file header
%   information for OBJ (the SEISMICFILEREADER object) according to the
%   number of traces in the scalar or vector TRACES. OFFSET is
%   the byte offset from the beginning of file to begin the read from.
%
%   Example:
%   % read in part of the trace header for a SEG2 formatted file.
%   s = SeismicFileReader('WFLT0001')
%
%   % read in the file header information to get the pointers to the trace
%   % data
%   fmt = {'uint16',{'FileDescriptor','RevisionNumber',...
%                     'TracePointerSize','NumberOfTraces'};
%          'uint8', {6, 'Terminators'} };
%   fileHeader = readFileHeader(s,fmt)
%
%   % get the number of traces and the trace locations in the file
%   fmt = {'uint32', {fileHeader.NumberOfTraces, 'TracePointers'}};
%   offset = 32; % bytes from beginning of file
%   locInFile = readFileHeader(s,fmt,offset);
%
%   % read in the header for trace 1
%   fmt = {'uint16',{'TraceDescriptor','TraceBlockSize'};...
%          'uint32',{'TraceDataSize','NumberOfSamples'};...
%          'uint8', {'TraceDataFormatCode'} };
%   offset = locInFile.TracePointers(1);      
%   traceHeader = readTraceHeader(s,fmt,offset)
%
%   % read in the data for trace 1
%   fmt = {'int32', traceHeader.NumberOfSamples}; % 32-bit integer 
%   offset = 1348; % trace 1 starts at byte 1348
%   traceData = readTraceData(s,fmt,offset)
%
%   See also SEISMICFILEREADER, READFILEHEADER, READTRACEHEADER

            % trace data should be a single numeric array
            if ~isOpen(obj)
                % open file if not already
                openFile(obj)
            end
            % fmt is assumed to be a 1x2 or 2x1 vector
            if ~isvector(fmt) || length(fmt) ~= 2
                error('SeismicFileReader:InvalideTraceFormat',...
                    'Format specification must be a vecor of length 2')
            end
            if nargin > 2
                offset = varargin{1};
                if offset >= 0
                    fseek(obj.FileID,offset,'bof');
                end % otherwise do nothing, use current position.
            end
            d = readFile(obj,fmt{2},fmt{1});
        end % readTraceData

        function st = closeFile(obj)
%CLOSEFILE closes the file
%   TF = CLOSEFILE(OBJ) closes the file that is linked to OBJ.  Returns 0
%   if close is successful, -1 if not.
            st = fclose(obj.FileID);
            obj.FileID = [];
        end
        
        function delete(obj)
%DELETE delete a SeismicFileReader object.
%   DELETE(OBJ) deletes the object OBJ, and closes any open files.

            % Make sure to close the file when deleted
            if ~isempty(obj.FileID)
                fclose(obj.FileID);
            end
        end % delete
        
    end % methods
    
    % Private Methods
    methods (Access = private)
        
        function openFile(obj,fileName)
        %OPENFILE open the file for reading
        %   OPENFILE(OBJ,FILENAME) opens the file in FILENAME and returns
        %   the file id to OBJ.FILEID.  An error is returned if the file
        %   cannot be opened.
            if nargin == 1
                fileName = fullfile(obj.FilePath,obj.FileName);
            end
            
            % If filename is provided, open it for reading if exists
            if ~isempty(fileName)
                assert(exist(fileName,'file') == 2, ...
                    'SeismicFileReader:FileNotFound',...
                    'Could not locate file: %s', fileName)
                
                obj.FileID = fopen(fileName,'r');
                
                assert(obj.FileID > 0, ...
                    'SeismicFileReader:CouldNotOpenFile',...
                    'Could not open file: %s', fileName);
            end
            
            % Set object properties from file, use which to capture path if
            % file is on the path
            
            ff = which(fileName);
            if isempty(ff)
                ff = fileName;
            end
               [obj.FilePath, obj.FileName, fileExt] = fileparts(ff);
               obj.FileName = [obj.FileName fileExt];
               s = dir(ff);
               obj.FileSize = s.bytes;
        end % openFile
        
        function st = isOpen(obj)
        %ISOPEN returns true or false if the file is open.
        %   ISOPEN(OBJ) returns TRUE if the file in OBJ.FILENAME is open,
        %   otherwise it returns FALSE.
        
            % Check if the file ID is valid
            try
                ferror(obj.FileID);
                st = true;
            catch e
                if sum(strcmpi(e.identifier,{'MATLAB:FileIO:InvalidFid',...
                        'MATLAB:badfid_mx'}))
                    st = false;
                else
                    throw(e)
                end
            end 
        end % isOpen
        
        function f = readFile(obj,sizeA,precision)
%READFILE reads format similar to fread
%   F = READFILE(OBJ,SIZEA,PRECISION) reads in data from the file linked to
%   OBJ using SIZEA and PRECISION specifiers (see FREAD for details).
%
%   See also FREAD
            if ~isOpen(obj)
                % open file if not already
                openFile(obj)
            end
            skip = 0;
            f = fread(obj.FileID,sizeA,precision,skip,obj.MachineFormat);
        end
        
        function h = readHeader(obj,fmt,varargin)
            if ~isOpen(obj)
                % open file if not already
                openFile(obj)
            end
            
            if nargin > 2
                offset = varargin{1};
                if offset >= 0
                    fseek(obj.FileID,offset,'bof');
                end % otherwise do nothing, use current position.
            elseif isempty(varargin)
                frewind(obj.FileID);
            end
            r = size(fmt,1);
            for i = 1:r
                nvars = length(fmt{i,2});
                % check if first index is a number
                if isnumeric(fmt{i,2}{1}) && nvars == 2
                    d = readFile(obj,fmt{i,2}{1},fmt{i});
                    h.(fmt{i,2}{2}) = d;
                else
                    d = readFile(obj,nvars,fmt{i});
                    for j = 1:nvars
                        h.(fmt{i,2}{j}) = d(j);
                    end % j loop
                end % if
            end % i loop
            
        end % readHeader
        
    end % private methods
    
end % class SeismicFileReader