classdef travelTimeMemmap < memmapfile   
    
    methods
        %% CONSTRUCTOR
        function obj = travelTimeMemmap(filename,sizeV,nRecords)
            obj = obj@memmapfile(filename,...
                'Format',{'double', [sizeV nRecords],'TT'},...
                'Repeat',1);
        end % travelTimeMemmap
        
        %% Sub-reference
        function varargout = subsref(obj, s)
            if strcmp(s(1).type,'()')
                if length(s) > 1
                    error('TRAVELTIMEMEMMAP:InvalidIndexing','This index operation is not allowed');
                end
                s2(1).type = '.';  s2(1).subs = 'Data';
                s2(2).type = '.';  s2(2).subs = 'TT';
                arrSize = get(obj,'format');
                arrSize = arrSize{2};
                offsetBytes = prod(arrSize);
                nRecords = arrSize(3);
                
                for i = 1:3
                    if strcmpi(s.subs{i},':');                        
                        arrSize(i) = arrSize(i);
                    else
                        arrSize(i) = length(s.subs{i});
                    end
                end
                data = zeros(arrSize);
        
                idx = 1;
                so = struct('type','.','subs','offset');
                for i = s.subs{3}
                    k = floor(i/nRecords);
                    if i > nRecords*k
                        % move offet
                        obj = subsasgn(obj,so,offsetBytes*8*k);
                        % adjust page index accordingly
                        p = i-k*nRecords;
                    elseif i == nRecords*k
                        % move offet
                        obj = subsasgn(obj,so,offsetBytes*8*(k-1));
                        % adjust page index accordingly
                        p = nRecords;
                    else
                        p = i;
                    end
                    
                    s2(3).type = '()'; s2(3).subs = [s.subs(1:2) p];
                    try
                        data(:,:,idx) = subsref@memmapfile(obj,s2);
                    catch me
                        if strcmpi(me.identifier,'MATLAB:memmapfile:unsupportedCSL')
                            % adjust format for too small of a file case
                            fmt = obj.Format;
                            tmp = fmt;
                            fmt{2}(3) = p;
                            obj.Format = fmt;
                            data(:,:,idx) = subsref@memmapfile(obj,s2);
                            obj.Format = tmp;
                        else
                            rethrow(me)
                        end % try-catch
                    end
                    %imagesc(data(:,:,idx)); drawnow;
                    idx = idx+1;
                end
                varargout{1} = data;
                
            else
                varargout{1} = subsref@memmapfile(obj,s);
            end
        end % Subsref
        
        %% Size method
        % Return the size of the virtual array
        function arrSize = size(obj,dim)
            arrSize = get(obj,'format');
            arrSize = arrSize{2};
            if exist('dim','var')
                arrSize = arrSize(dim);
            end
        end % size
            
    
    end % Methods
    
end % travelTimeMemmap