classdef SegYFileReader < SeismicFileReader
%SegYFileReader SEG Y file reader class.
%
%   SegYFileReader is a SeismicFileReader class that provides basic
%   functionality for reading a SEG Y formatted file.  This class can read
%   SEG Y Version 0/1 files.
%
%   SegYFileReader properties:
%       FileDescriptor      - SEG2 file identifier label (must = 3A55h)
%       FileName            - name of the seismic binary file
%       FilePath            - path of the file
%       FileFormat          - file format description
%       FileHeader          - Header information from the seismic file
%       MachineFormat       - Binary file format, big/little endian
%       NumberOfTraces      - Number of traces in the file
%       TraceDescription    - Trace information from the seismic file
%
%   See also SEISMICFILEREADER, SEG2FILEREADER
    
    properties (SetAccess = protected)
        TapeLabel = '';

    end % public properties
    
    methods
        function obj = SegYFileReader(fileName,procHeader,procTraceHeader,endian)
        % Constructor for Seg2FileReader
            if ~exist('fileName','var')
                % creat an empty clase
                fileName = '';
                procHeader = false;
            end
            
            if ~exist('procHeader','var')
                % read in the header
                procHeader = true;
            end
            
            if ~exist('procTraceHeader','var')
                % read in the header
                procTraceHeader = true;
            end
            
            if ~exist('endian','var')
                endian = 'ieee-be'; % SEGY Rev 1.0 should be big endian
            end
            
            obj = obj@SeismicFileReader(fileName);
            obj.MachineFormat = endian;            
            
            if procHeader
                readFileHeader(obj);
            end
            
            if procHeader && procTraceHeader
                readTraceHeader(obj);
            end
        end
        
        function fh = readFileHeader(obj)
%READFILEHEADER reads SEG Y header information
%   H = READFILEHEADER(OBJ) reads in the optional tapel lablel header, the
%   3200 bytes textual header, the 400 byte binary header, and option 3200
%   byte textual headers.
%
%   Example: TODO
%
%   See also SEGYFILEREADER, READTAPEHEADER, READTAPEDATA

            % read in the optional header if present
            % Is there a optional file header?
            fmt = {'uchar',{128,'TapeLabel'}};
            fh = readFileHeader@SeismicFileReader(obj,fmt);
            tapeHeader = char(fh.TapeLabel');
            
            % if SYX.X is present, this file contains an optional header
            if strcmpi('SY',tapeHeader(5:6))
                obj.TapeLabel.StorageUnitSeqNumber = tapeHeader(1:4);
                obj.TapeLabel.SEGYRevision = tapeHeader(5:9);
                obj.TapeLabel.StorageUnitStructure = tapeHeader(10:15);
                obj.TapeLabel.BindingEdition = tapeHeader(16:19);
                obj.TapeLabel.MaxBlockSize = tapeHeader(20:29);
                obj.TapeLabel.ProducerOrganizationCode = tapeHeader(30:39);
                obj.TapeLabel.CreationDate = tapeHeader(40:50);
                obj.TapeLabel.SerialNumber = tapeHeader(51:62);
                obj.TapeLabel.Reserved = tapeHeader(63:68);
                obj.TapeLabel.StorageSetIdentifier = tapeHeader(69:128);
            else
                obj.TapeLabel = 'None';
            end
            
            % if present, set the offset to account for it
            if isstruct(obj.TapeLabel)
                offset = 128;
            else
                offset = 0;
            end
            
            % read in the 3200 byte required header
            fmt = {'uchar',{3200,'TextHeader'}};
            fh = readFileHeader@SeismicFileReader(obj,fmt,offset);
            obj.FileHeader.TextualHeader = char(fh.TextHeader');
            
            % read in the 400 Byte binary header
            offset = offset + 3200;
            fmt = {'int32',{'JobID','LineNumber','ReelNumber'};
                   'int16',{'DataTracePerEnsemble',...
                            'AuxillaryTracePerEnsemble',...
                            'SamplingTimeInMicroSec',...
                            'OrigninalSamplingTime',...
                            'NumberOfSamples','OriginalNumberOfSamples',...
                            'SampleFormatCode','EnsembleFold',...
                            'TraceSorting','VerticalSumCode',...
                            'SweepFreqencyStart','SweepFrequencyEnd',...
                            'SweepLength','SweepType','SweepChannel',...
                            'SweepTaperLengthStart',...
                            'SweepTaperLengthEnd','TaperType',...
                            'CorrelatedTraces','BinaryGain',...
                            'AmplitudeRecovery','MeasurementSystem',...
                            'ImpulsePolarity','VibratoryPolarity'};
                    'int16',{120,'Unassigned1'};
                    'int16',{'RevisionNumber','FixedLengthTraceFlag',...
                             'NumberOfExtendedTextualHeaders'};
                    'int16',{94,'Unassigned2'}
                    };
                
             fh = readFileHeader@SeismicFileReader(obj,fmt,offset);
             obj.FileHeader.BinaryHeader = fh;
             
             % Verify SEG Y supported revision
             v = obj.FileHeader.BinaryHeader.RevisionNumber;
             if v > 257
                 error('SegYFileReader:UnsupportedVersion',...
                     'Version number found in file is %2.2f. Supported versions are up to 1.0',...
                     num2str(v/10));
             else
                 obj.FileFormat = ['SEG Y Version ',num2str(v/10)];
             end
             
             % Set number of traces
             obj.NumberOfTraces = obj.FileHeader.BinaryHeader.DataTracePerEnsemble;
                             
             offset = offset+400;
             
             % Read Extended 3200 byte textual headers if present
             for i = 1:obj.FileHeader.BinaryHeader.NumberOfExtendedTextualHeaders;
                 fmt = {'uchar',{3200,'TextHeader'}};
                 fh = readFileHeader@SeismicFileReader(obj,fmt,offset);
                 obj.FileHeader.ExtendedTextualHeader(i) = char(fh.TextHeader');
                 offset =offset+3200;
             end
             
             
             % Validate input selections
             switch obj.FileHeader.BinaryHeader.SampleFormatCode
                 case 1 % 4-byte IBM floating-point
                     obj.FileHeader.TraceFormat = 'uint';
                     mult = 4;
                 case 2 % 4-byte two's complement integer
                     obj.FileHeader.TraceFormat = 'int32';
                     mult = 4;
                 case 3 % 2-byte two's complement integer
                     obj.FileHeader.TraceFormat = 'int16';
                     mult = 2;
                 case 4 % 4-byte fixed point with gain (obsolete)
                     obj.FileHeader.TraceFormat = 'int32';
                     mult = 4;
                 case 5 % 4-byte IEEE floating-point
                     obj.FileHeader.TraceFormat = 'float32';
                     mult = 4;
                 case 8 % 1-byte two's complement integer
                     obj.FileHeader.TraceFormat = 'int8';
                     mult = 1;
                 otherwise
                     warning('Couldn''t determine data format, assuming float 32');
                     obj.FileHeader.TraceFormat = 'float32';
                     mult = 4;
             end
             
             % Store trace offset pointers
             fmt = {'int16',{1,'NumberOfSamples'}};
             fh = readFileHeader@SeismicFileReader(obj,fmt,offset+114);
             traceSize = fh.NumberOfSamples*mult + 240;
             mxTraces = ceil((obj.FileSize - offset)/traceSize);
             obj.FileHeader.TracePointers = zeros(mxTraces,1);
             obj.FileHeader.TracePointers(1) = offset;
             % check to see if all traces are of equal size
             if obj.FileHeader.BinaryHeader.FixedLengthTraceFlag == 1
                 % all are the same size so build pointers
                 i = 1:(mxTraces+1);
                 obj.FileHeader.TracePointers = [(i-1)*traceSize + offset]';
                 obj.NumberOfTraces = mxTraces;
             else % traces could be different sizes
                i = 1;
                 while (offset<obj.FileSize) %for i = 2:obj.NumberOfTraces;
                     % read the trace headers to get the size of the traces
                     % number of samples is 114 bytes into the trace header
                     i = i+1;
                     fh = readFileHeader@SeismicFileReader(obj,fmt,offset+114);
                     ptr = fh.NumberOfSamples*mult + 240 +...
                            obj.FileHeader.TracePointers(i-1);
                     if ptr < obj.FileSize
                         obj.FileHeader.TracePointers(i) = ptr;
                     end
                    offset = ptr;
                 end
                 obj.NumberOfTraces = i;
                 if i < mxTraces
                    obj.FileHeader.TracePointers(i+1:end) = [];
                 end
             end

        end %readFileHeader
        
        function th = readTraceHeader(obj,traceNumber)
            if nargin < 2
                traceNumber = 1:obj.NumberOfTraces;
            end
            % read in header data
            fmt = {'int32',{'TraceSequenceLine','TraceSequenceFile',...
                            'OriginalFieldRecord','TraceNumberInField',...
                            'EnergySourceNumber','EnsembleNumber',...
                            'EnsembleTraceNumber'};
                   'int16',{'TraceIdentificationCode','NumberVertSummedTraces',...
                            'NumberHorizSummedTraces','DataUse'};
                   'int32',{'ReceiverGroupOffset',...
                            'ReceiverGroupElevation',...
                            'SurfaceSourceElevation',...
                            'SurfaceSourceDepth',...
                            'ReceiverDatumElevation',...
                            'SourceDatumElevation',...
                            'SourceWaterDepth','GroupWaterDepth'};
                    'int16',{'ElevationOrDepthScalar','CoordinateScalar'};
                    'int32',{'SourceX','SourceY','GroupX','GroupY'};
                    'int16',{'CoordinateUnits','WeatheringVelocity',...
                             'SubweatheringVelocity',...
                             'SourceUpholeTime_ms','GroupUpholeTime_ms',...
                             'SourceStaticCorrection_ms',...
                             'GroupStaticCorrection_ms',...
                             'TotalStaticApplied_ms',...
                             'ALagTime_ms','BLagTime_ms',...
                             'RecordingDelay_ms',...
                             'MuteTimeStart_ms','MuteTimeEnd_ms',...
                             'NumberOfSamples',...
                             'SampleInterval_ms',...
                             'GainType','GainConstant_dB',...
                             'InitialGain_dB','Correlated',...
                             'SweepFrequencyStart_Hz',...
                             'SweepFrequencyEnd_Hz',...
                             'SweepLength_ms','SweepType',...
                             'SweepTaperLengthStart_ms',...
                             'SweepTaperLengthEnd_ms',...
                             'TaperType','AliasFilter_Hz',...
                             'AliasSlope_dBperOctave',...
                             'NotchFilterFrequency_Hz',...
                             'NotchFilterSlope_dBperOctave',...
                             'LowCutFrequency_Hz',...
                             'HighCutFrequency_Hz',...
                             'LowCutSlope_dBperOctave',...
                             'HighCutSlope_dBperOctave',...
                             'Year','Day','Hour','Minute','Second',...
                             'TimeBasisCode','WeightingFactor',...
                             'GeophoneGroupNumberRoll1',...
                             'GeophoneGroupNumberFirstTraceOriginalRecord',...
                             'GeophoneGroupNumberLastTraceOriginalRecord',...
                             'GapSize','OverTravel'};
                    'int32',{'EnsembleX','EnsembleY','InLineNumber3D',...
                             'CrossLineNumber3D','ShotPointNumber'};
                    'int16',{'ShotPointScalar','TraceMeasurementUnit'};
                    'int32',{1,'TransductionConstantMantissa'};
                    'int16',{'TransductionConstantExponent',...
                             'TranductionUnits','TraceIdentifier',...
                             'TimeScalar','SourceTypeOrientation'};
                    'bit48',{1,'SourceEnergyDirection'};
                    'int32',{1,'SourceMeasurementMantissa'};
                    'int16',{'SourceMeasurementExponent',...
                             'SourceMeasurementUnit'};
                    'int16',{4,'Unassigned'}
                    };
                       
            for t = 1:length(traceNumber)
                offset = obj.FileHeader.TracePointers(traceNumber(t));      
                f = readTraceHeader@SeismicFileReader(obj,fmt,offset);
                th(t) = f;
            end
            
            if nargout == 0 % add to object if no output requested
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
                
                %if  th(t).NumberOfSamples > r
                %    th(end:end+n,:) = NaN;
                %end 
                    
                % get trace data and fill pre-allocated array
                fmtCode = obj.FileHeader.TraceFormat;
                fmt = {fmtCode, th(1).NumberOfSamples};
                offset = obj.FileHeader.TracePointers(traceNumber(t)) + ...
                          240; % offset for trace binary header
                % account for ibm floating point
                if obj.FileHeader.BinaryHeader.SampleFormatCode == 1
                    td(:,t) = ibm2ieee(readTraceData@SeismicFileReader(obj,fmt,offset));
                else
                    td(:,t) = readTraceData@SeismicFileReader(obj,fmt,offset);
                end
            end % for loop
        end % readTraceData
        
        function txt = parseTextBlock(obj,textBlock)
            % parse out text block data
            txt = [];
                % extract keywords 
                [keys,idx] = regexp(textBlock, '[A-Z_]\w*', 'match','start');
                % check validity
                [~,idx] = obj.validKeys(keys,idx);
                
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
   
        function varargout = subsref(obj,s)
            if numel(obj) > 1
                if length(s) > 1
                    obj = obj(s(1).subs{:});
                    varargout{1} = subsref(obj,s(2:end));
                else
                    varargout{1} = obj(s.subs{:});
                end
            else
                if strcmp(s(1).type,'()')
                    if length(s) > 1
                        error('SeismicFileReader:InvalidIndexing','This index operation is not allowed');
                    end
                    if strcmpi(s.subs{2},':')
                        ind = [];
                    else
                        ind = s.subs{2};
                    end
                    tmp = readTraceData(obj,ind);
                    varargout{1} = tmp(s.subs{1},:);
                elseif strcmp(s.type,'.')
                    varargout{1} = obj.(s.subs);
                end
            end
        end % subsref
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
