classdef TraceData < File
    % Class for SEGY Trace Headers
    % Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
       data;
    end
    
    methods 
         function obj = TraceData(file, varargin)
             obj=obj@File(file,varargin{:});
         end
         
         
    end
    
    methods (Static)
    end
end



