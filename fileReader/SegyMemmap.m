classdef SegyMemmap < memmapfile
    
    methods
        %% CONSTRUCTOR
        function obj = SegyMemmap(filename,nSamples,nTraces)
            obj = obj@memmapfile(filename,...
                'Offset',3600,...
                'Format',{'int32'  [60 1] 'traceHeader';...
                'single' [nSamples 1] 'trace'},...
                'Repeat',nTraces);
        end % SegyMemmap
        
        %% Sub-reference
        function varargout = subsref(obj, s)
            if s(1).type == '()'
                if length(s) > 1
                    error('SEGYMEMAP:InvalidIndexing','This index operation is not allowed');
                end
                mxc = length(s.subs{2});
                mxr = length(s.subs{1});
                getAllCol = false;
                if strcmpi(s.subs(2),':')
                    s2(1).type = '.'; s2(1).subs = 'Repeat';
                    mxc = subsref@memmapfile(obj,s2);                    
                    getAllCol = true;
                end
                if strcmpi(s.subs(1),':')
                    s2(1).type = '.'; s2(1).subs = 'Data';
                    s2(2).type = '()'; s2(2).subs = {1};
                    s2(3).type = '.'; s2(3).subs = 'trace';
                    mxr = length(subsref@memmapfile(obj,s2));
                end
                
                s2(1).type = '.';  s2(1).subs = 'Data';
                s2(3).type = '.';  s2(3).subs = 'trace';
                s2(4).type = '()'; s2(4).subs = s.subs(1);
                data = zeros(mxr,mxc);
                for c = 1:mxc
                    s2(2).type = '()';
                    if getAllCol
                        s2(2).subs = {c};
                    else
                        s2(2).subs = {s.subs{2}(c)};
                    end
                    data(:,c) = swapbytes(subsref@memmapfile(obj,s2));
                end
                varargout{1} = data;
            else
                data = subsref@memmapfile(obj,s);
                if strcmpi(s(1).subs,'data') && ~isa(data,'struct')
                    varargout{1} = swapbytes(data);
                else
                    varargout{1} = data;
                end
            end
        end % Subsref
    end % Methods
end % SegyMemmap