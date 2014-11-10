classdef HeaderDefinitions < File
    %Class for SEGY header definitions
    properties
        keys;
        values;
        % column indices
        n; %name
        s; %startbyte
        e; %endbyte
        t; %type
        r; %recommended
        d; %description
        c; %comment
        i; %interpretation
    end

    methods
        function obj = HeaderDefinitions(file, varargin)
            obj=obj@File(file,varargin{:});
            if (nargin>0)
                obj.readDefinitions;
                %get column index for header definitions.
                obj.n = strmatch(upper('name'),upper(obj.keys));
                obj.s = strmatch(upper('startbyte'),upper(obj.keys));
                obj.e = strmatch(upper('endbyte'),upper(obj.keys));
                obj.t = strmatch(upper('type'),upper(obj.keys));
                obj.r = strmatch(upper('recommended'),upper(obj.keys));
                obj.d = strmatch(upper('description'),upper(obj.keys));
                obj.c = strmatch(upper('comment'),upper(obj.keys));
                obj.i = strmatch(upper('interpretation'),upper(obj.keys));
            end
        end
    end

end

