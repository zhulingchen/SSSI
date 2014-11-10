classdef Trace 
    % Class for SEGY Trace Headers
    % Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        traceheader;
       tracedata;
    end
    
    methods 
         function obj = Trace(obj_file, varargin)
             if ischar(obj_file)
             obj.traceheader=TraceHeader(obj_file,varargin{:});
            obj.tracedata=TraceData(obj_file,varargin{:});
             elseif isa(obj_file,'TraceHeader')
                 obj.traceheader=obj_file;
                 obj.tracedata=TraceData(obj.traceheader.filename,varargin{:});
             end
         end
         function obj=getTraces(obj,varargin)
             if isempty(obj.traceheader.traceoffsets) && isempty(obj.tracedata.data)
                 obj=obj.getTracesfromfile(varargin{:});
             else
         obj=obj.getTraceswithtracehead(varargin{:});
             end
         end
         function writetrace(obj)
            obj.writeTrace();
         end
         function delete(obj)
             obj = obj.closeFile();            
         end
    end
    
    methods (Static)
        [type bytes endian]=uigettracetype(obj,numoftraces);
    end
end



