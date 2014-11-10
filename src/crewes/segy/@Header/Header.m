classdef Header < File
    % Class for SEGY Trace Headers
    % Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        header;
        definitions;
        hdrsize; %bytes
        hdroffset; %bytes from beginning of file.
        nontypecasthdr=uint8([]); % raw uint8 header data
        filefmt='';
    end
    
    methods
        function obj = Header(file, varargin)
            obj=obj@File(file,varargin{:});
        end
        function obj.setheadervalue(obj,word,value,scalevalues)
            obj=obj.setheadervalue(obj,word,value,scalevalues);
        end
        function obj=obj.uiHeaderConverter(obj)
            obj=obj.uiHeaderConverter(obj);
        end
        function obj=obj.standardizeheader(obj)
            obj=standardizeheader(obj);
        end
    end
    
    methods (Static)
        val=obj.getheadervalue(obj,word,scalevalues);
        endi=getfileendiantype4sgy(obj,offset);
        traceheaderposition=uitestforextendedheaders(obj);
        %[values colnames] = readHeaderDefinitions(filename);
        %this needs much work
    end
end





