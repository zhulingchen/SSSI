classdef Seg2FileReader < SeismicFileReader
%Seg2FileReader SEG2 file reader class.
%
%   Seg2FileReader is a SeismicFileReader class that provides basic
%   functionality for reading a SEG2 formatted file.  This class can read
%   SEG2 Version 0 files.
%
%   Seg2FileReader properties:
%       FileDescriptor      - SEG2 file identifier label (must = 3A55h)
%       FileName            - name of the seismic binary file
%       FilePath            - path of the file
%       FileFormat          - file format description
%       FileHeader          - Header information from the seismic file
%       MachineFormat       - Binary file format, big/little endian
%       NumberOfTraces      - Number of traces in the file
%       TraceDescription    - Trace information from the seismic file
%
%   See also SEISMICFILEREADER, SEGYFILEREADER
    
    properties (SetAccess = protected)
        
%FILEDESCRIPTOR Identifier label.
%   The FILEDESCRIPTOR property records the file identifier.  For SEG2
%   formatted files, this must be hexidecimal 3A55.  This number is used to
%   identify the file as SEG2 and determine if the file is using big or
%   little ENDIAN format.

        FileDescriptor = '';

    end % public properties
    
    methods
        function obj = Seg2FileReader(fileName,readHeader)
            % Constructor for Seg2FileReader
            if nargin == 0 % create an empty class
                fileName = '';
                readHeader = false;
            elseif nargin == 1 
                % read file header unless told not to
                readHeader = true;
            end
            
            obj = obj@SeismicFileReader(fileName);
            
            if readHeader
                readFileHeader(obj);
                readTraceHeader(obj);
            end
        end
        
        function fh = readFileHeader(obj)
            % file header reader
            fmt = {'uint8',{'byte1','byte2'}};
            % Determine if this is a Seg2 formatted file and machine format
            % type
            fh = readFileHeader@SeismicFileReader(obj,fmt);
            if fh.byte1 < fh.byte2
                fileDescriptor = [dec2hex(fh.byte1),dec2hex(fh.byte2),'h'];
            else
                fileDescriptor = [dec2hex(fh.byte2),dec2hex(fh.byte1),'h'];
            end
            e = MException('Seg2FileReader:InvalidFormat',...
                'Unrecognized SEG2 format\n\t First byte in record is %sh\n\t SEG2 (1990) requires 3Ah (big-endian) or 55h (little-endian)\n\t File descriptor block value: %s',...
                dec2hex(fh.byte1),fileDescriptor);
            switch upper(dec2hex(fh.byte1))
                case '3A'
                    obj.MachineFormat = 'ieee-be'; % big-endian
                case '55'
                    obj.MachineFormat = 'ieee-le'; % little-endian
                otherwise
                    throw(e)
            end % switch
            
            % double check the 1st two bytes are indeed 3A55h
            if ~strcmpi(fileDescriptor,'3A55h')
                throw(e)
            end
            obj.FileDescriptor = fileDescriptor;
            
            % read in header data 
            fmt = {'uint16',{'FileDescriptor','RevisionNumber',...
                             'TracePointerSize','NumberOfTraces'};
                   'uint8', {'StringTerminatorSize',...
                             'FirstStringTerminator',...
                             'SecondStringTerminator',...
                             'SizeofLineTerminator',...
                             'FirstLineTerminator',...
                             'SecondLineTerminator'};
                   };
            fh = readFileHeader@SeismicFileReader(obj,fmt);  
      
            % set the file format
            obj.FileFormat = ['SEG2 Revision ',num2str(fh.RevisionNumber),...
                              ' ',obj.MachineFormat];
                          
            % set the number of traces
            obj.NumberOfTraces = fh.NumberOfTraces;
            
            % read in variable length header data
            fmt = {'uint32', {obj.NumberOfTraces, 'TracePointers'}};
            offset = 32; % bytes from beginning of file
            tmp = readFileHeader@SeismicFileReader(obj,fmt,offset);
            fh.TracePointers = tmp.TracePointers;
            
            % read in text block
            fmt = {'uint8', {fh.TracePointers(1)-32-fh.TracePointerSize,...
                             'TextBlock'}};
            offset = 32+fh.TracePointerSize; % bytes from beginning of file
            tmp = readFileHeader@SeismicFileReader(obj,fmt,offset);
                      
            % remove nonprintable characters
            fh.RawTextBlock = char(tmp.TextBlock');
            fh.TextData = obj.parseTextBlock(fh.RawTextBlock);
            obj.FileHeader = fh;
            
        end %readFileHeader
        
        function th = readTraceHeader(obj,traceNumber)
            if nargin < 2
                traceNumber = 1:obj.NumberOfTraces;
            end
            % read in header data
            fmt = {'uint16',{'TraceDescriptor','TraceBlockSize'};...
                   'uint32',{'TraceDataSize','NumberOfSamples'};...
                    'uint8', {'TraceDataFormatCode'} };
         
            for t = 1:length(traceNumber)
                offset = obj.FileHeader.TracePointers(traceNumber(t));      
                f = readTraceHeader@SeismicFileReader(obj,fmt,offset);
                
                % check to make sure we are reading a header
                assert(strcmpi(dec2hex(f.TraceDescriptor),'4422'),...
                      'Seg2FileReader:InvalidTraceFormat',...
                      'Unrecognized SEG2 format\n\t First byte in trace %d record is %sh\n\t SEG2 (1990) requires 4422h\n\t',...
                      t,dec2hex(f.TraceDescriptor));
                  
                % read in text block            
                fmt2 = {'uint8', {f.TraceBlockSize - 32, 'TextBlock'}};
                offset = 32+offset; % bytes from beginning of file
                tmp = readTraceHeader@SeismicFileReader(obj,fmt2,offset);
                % remove nonprintable characters
                %tmp.TextBlock(tmp.TextBlock < 32 & tmp.TextBlock > 0) = [];
                f.RawTextBlock = char(tmp.TextBlock');
                f.TextData = obj.parseTextBlock(f.RawTextBlock);
                th(t) = f;
            end
            
            if nargout == 0
                obj.TraceDescription = th;
            end
        end % readTraceHeader
        
        function [td,th] = readTraceData(obj,traceNumber,td)
            if nargin == 1 || isempty(traceNumber)
                traceNumber = 1:obj.NumberOfTraces;
            end
            
            % get 1st trace size to determine initial array size
            th = readTraceHeader(obj,1);
            r = th.NumberOfSamples;
            c = length(traceNumber);
            
            % create as distributed array if requested
            if nargin > 2 && license('test','distrib_computing_toolbox')
                td = distributed.nan(r,c);
            else
                try
                    td = nan(r,c);
                catch
                    e = lasterror;
                    if strcmpi(e.identifier,'MATLAB:nomem')
                        % test maximum array size
                        u = memory;
                        str = ['\n\tCould not allocate an array of ',...
                            num2str(r*c*8),'bytes.  Available memory is: \n',...
                            u.MaxPossibleArrayBytes,...
                            '\tTry declaring as a distributed array:\n\t\t',...
                            'd = readTraceData(obj,traceNumbers,distributed)'];
                        e.message = sprintf([e.message,str]);
                    end
                   rethrow(e)
                end
                        
                    
            end
           
            for t = 1:length(traceNumber)
                % read the trace header information
                th(t) = readTraceHeader(obj,traceNumber(t));
                
                % need to account for variable length data (NaNs)
                n = th(t).NumberOfSamples;
                % since r is sized by the first data set, it should never
                % be smaller than this, we only need to test for the case
                % where it is greater
                
                if  th(t).NumberOfSamples > r
                    th(end:end+n,:) = NaN;
                end 
                    
                % get trace data and fill pre-allocated array
                switch th(t).TraceDataFormatCode
                    case 1
                        fmtCode = 'int16'; % 16-bit fixed point
                    case 2
                        fmtCode = 'int32'; % 32-bit fixed point   
                    case 3
                        fmtCode = 'real*20'; % 20-bit floating point
                    case 4
                        fmtCode = 'float32'; % 32-bit floating point
                    case 5
                        fmtCode = 'float64'; % 64-bit floating point
                    otherwise
                        e = MException('Seg2FileReader:InvalidDataFormat',...
                            'Unrecognized SEG2 format\n\t Trace data format code in record %d is %d\n\t SEG2 (1990) specifies 1-5 format code\n\t',...
                            t,th(t).TraceDataFormatCode);
                        throw(e);
                end % switch
                fmt = {fmtCode, th(1).NumberOfSamples};
                offset = obj.FileHeader.TracePointers(traceNumber(t)) + ...
                          th(t).TraceBlockSize;
                td(:,t) = readTraceData@SeismicFileReader(obj,fmt,offset);
            end % for loop
        end % readTraceData
        
        function txt = parseTextBlock(obj,textBlock)
            % parse out text block data
            txt = [];
                % extract keywords 
                [keys,idx] = regexp(textBlock, '[A-Z_]\w*', 'match','start');
                % check validity
                [tmp,idx] = obj.validKeys(keys,idx);
                
                % remove 1st two bytes (offsets) in file from keywords
                % and last offset bytes
                textBlock([idx-2, idx-1,end-1:end]) = [];
                                
                % extract keywords with new index locations 
                [keys,idx] = regexp(textBlock, '[A-Z_]\w*', 'match','end');
                % check validity
                [keys,idx] = obj.validKeys(keys,idx);

                for i = 1:length(idx)
                    if i < length(idx)
                        str = textBlock((idx(i)+1):(idx(i+1)-length(keys{i+1})));
                    else
                        str = textBlock((idx(i)+1):end);
                    end
                    % check if it is a date or time string
                    if isempty(strfind(str,':')) && isempty(strfind(str,'\'))
                        tmpStr = str2num(str);
                    else
                        tmpStr = deblank(strtrim(str));
                    end
                    if isempty(tmpStr)
                        % If only whitespaces don't add
                        if ~isempty(strtok(str,char(0)))
                            txt.(keys{i}) = deblank(strtrim(str));
                        end
                    else
                        txt.(keys{i}) = tmpStr;
                    end
                end % for
        end % parseTextBlock
   
    end % methods
    
    methods (Static)
        
        function [keys,idx] = validKeys(keys,idx)
            % check valid keys
            validKey = {'ACQUISITION_DATE',... % file header
                        'ACQUISITION_TIME',...
                        'CLIENT',...
                        'COMPANY',...
                        'GENERAL_CONSTANT',...
                        'INSTRUMENT',...
                        'JOB_ID',...
                        'OBSERVER',...
                        'PROCESSING_TIME',...
                        'TRACE_SORT',...
                        'UNITS',...
                        'NOTE',...
                        'ALIAS_FILTER',... % trace header
                        'AMPLITUDE_RECOVERY',...
                        'BAND_REJECT_FILTER',...
                        'CDP_NUMBER',...
                        'CDP_TRACE',...
                        'CHANNEL_NUMBER',...
                        'DATUM',...                       
                        'DELAY',...
                        'DESCALING_FACTOR',...
                        'DIGITAL_BAND_REJECT_FILTER',...
                        'DIGITAL_HIGHT_CUT_FILTER',...
                        'DIGITAL_LOW_CUT_FILTER',...
                        'END_OF_GROUP',...
                        'FIXED_GAIN',...
                        'HIGH_CUT_FILTER',...
                        'LINE_ID'...
                        'LOW_CUT_FILTER',...
                        'NOTCH_FREQUENCY',...
                        'POLARITY',...
                        'RAW_RECORD',...
                        'RECEIVER',...
                        'RECEIVER_GEOMETERY',...
                        'RECEIVER_LOCATION',...
                        'RECEIVER_SPECS',...
                        'RECEIVER_STATION_NUMBER',...
                        'SAMPLE_INTERVAL',...
                        'SHOT_SEQUENCE_NUMBER',...
                        'SKEW',...
                        'SOURCE',...
                        'SOURCE_GEOMETRY',...
                        'SOURCE_LOCATION',...
                        'SOURCE_STATION_NUMBER',...
                        'STACK',...
                        'STATIC_CORRECTIONS',...
                        'TRACE_TYPE',...
                        'UNIT',...% non-spec keywords below                   
                        'ACQUISITION_SECOND_FRACTION'};
            
             lidx = ismember(keys,validKey);
             keys = keys(lidx);
             idx = idx(lidx);
        end % validKeys
                      
    end % Static methods
    
end % Seg2FileReader Class
