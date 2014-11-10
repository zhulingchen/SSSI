classdef BinaryHeader < Header
    % Class for SEGY Binary Header
    %   BinaryHeader contains the binary header for the SEG-Y file.  This
    %     header is stored as uint8 type values.  The definition file
    %     directs matlab how to interpret these values therefor it is
    %     essential that getheadervalue and setheadervalue are used.
    %
    %  When BinaryHeader is called filename must be entered to construct a
    %    new BinaryHeader object, other arguments can be entered as listed
    %    below.
    %
    % Arguments that can be entered into traceheader are:
    %   'bindefinition', this will allow the binaryheader to be decoded using
    %      a non-standard definitions file.  'bindefinition' must be
    %      accompanied by a *.csv filename
    %   'new', this is a flag that will create a traceheader object without
    %      reading in any of the file.  This is especially usefull for when 
    %      a new file is created. 'new' is accompanied by '1'.
    %
    %

    properties
    
    end

    methods
        function obj = BinaryHeader(file, varargin)
            obj=obj@Header(file,varargin{:});
            
            obj.hdrsize=400; %bytes
            obj.hdroffset=3200; %bytes from beginning of file.
            guessByteOrder(obj);
            new=0;
            definitions='segyrev1_file.csv';
            for i = 1:2:length(varargin)
                name = varargin{i};
                value = varargin{i+1};
                if strcmpi(name,'machineformat');
                    obj.machineformat=value;
                elseif strcmpi(name,'new')
                    new=1; 
                    elseif strcmpi(name,'bindefinitions')
                    definitions=value;  
                end
            end
            if (nargin>0)
                obj.definitions = HeaderDefinitions(definitions);
                info=dir(obj.filename);
                if info.bytes<=obj.hdroffset || new
                    obj.nontypecasthdr=uint8(zeros(obj.hdrsize,1));
                    [ab, ac, e] = computer;
                    obj.filefmt=e;
                else
                obj.filefmt=obj.getfileendiantype4sgy(obj,3224);
                readHeader(obj);
                %convertHeader2(obj)
                
                end
            end
        end
        function obj=convertHeader(obj)
            obj=convertHeader2(obj);
        end
        function writeheader(obj)
            obj.writeHeader();
        end
        
        function delete(obj)
             obj.closeFile();            
         end
%         function bhead = readBinaryHeader(fileID);
%         writeBinaryHeader(obj);
%         getBinaryHeaderValue(obj);
%         setBinaryHeaderValue(obj);
%         interpretBinaryHeader(obj);
    end
end



