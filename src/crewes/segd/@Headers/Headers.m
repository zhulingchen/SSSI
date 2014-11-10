classdef Headers
    %Class for SEGD header definitions
    properties
        headerType;   %One of:
                      % generalHeader1
                      % generalHeader2
                      % generalHeaderN
                      % scanTypeHeader
                      % demuxTraceHeader
                      % traceHeaderExtension
                      % generalTrailer
        fid;          % Open file id
        filePosition; % bytes from start of file (for debugging)
        header;       % SEG-D header from disk as 8-bit unsigned integers
        hexHeader;    % SEG-D header converted to hex string
        headerSize;   % 32 bytes, except for trace headers
        headerDefCols;% header word definition column names
        headerDefRows;% header word definition row names
        headerDefs;   % header word definitions
    end

    methods
        function obj = Headers(fid, headerType, varargin)
            obj.fid = fid;
            obj.headerType = headerType;
            obj.headerSize = 32;
            obj.filePosition = ftell(obj.fid);
            if(strcmp(headerType,'demuxTraceHeader'))
                obj.headerSize = 20;
            end

            obj = obj.readHeaderDefinitions();
            obj = obj.readHeader();
        end
    end

end

