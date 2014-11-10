function saveascii(filename,varargin)
%SAVEASCII  Saves the 2-D array x into an ascii file, just like the
%built-in Matlab save command.  The builtin save doesn't work with the
%'-ascii' in compiled code.  This provides the same function. 
%
%Unlike the matlab save function, matrices to dump should be passed as
%data, not as strings:
%Example:
% a=[1 2 ; 3 4];
% save('file.txt','a','-ascii');     % Matlab built-in 
% saveascii('file.txt',a);           % this is equivalent
%
% b=[5 6 7 8];
% saveascii('file.txt','a','b','-ascii');
% saveascii('file.txt',a,b);         % same as above
%
% Author: Henry Bland  Nov 2002

    fp = fopen(filename,'w');
    if (fp <= 0) 
        error(['Unable to write to ' filename ]);
    end
    for arg=1:length(varargin)
        for r = 1:size(varargin{arg},1);
            x = varargin{arg};
            fprintf(fp,'%16.7e',x(r,:));
            fprintf(fp,'\n');
        end
    end 
    fclose(fp);  
    
            
