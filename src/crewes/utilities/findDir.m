function result = findDir(dirname,startdir,result)
%
% function result = finddir(dirname,startdir)
%
% Recursively looks for a directory name
%  Where: dirname  is the directoryname to be found (exact match)
%         startdir is the top level directory to start searching in        
%         result   is the result of the search
%  Usage: result = findDir(dirname)
%         result = findDir(dirname,startdir)
%

% Check for first call to findDir
if nargin <1 || nargin >3
   disp('Use ''help findDir'' for usage');
   return
elseif nargin == 1
    startdir = userpath();
    startdir = startdir(1:end-1);
    result=[];
elseif nargin == 2
    result=[];
end

%disp(startdir)

fl     = dir(startdir);   %get list of all files/dirs in userdir
n      = [fl.isdir];     %get logical indices for just the dirs
dirs   = { fl(n).name }; %cell array list of dir names


s = strcmp(dirs,dirname); %check for exact match
if any(s) %success!
    result = fullfile(startdir,dirs{s});
    return
else %failure, recurse
    for k = 1:length(dirs)
        if ~strncmp(dirs{k},'.',1) %skip directories '.' and '..'
            result = findDir(dirname,fullfile(startdir,dirs{k}),result);
            if ~isempty(result) %success!!!
                return
            end
        end
    end
end

end %end findDir()
