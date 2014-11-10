function val = toInterpretedString ( obj )
%
%function val = toInterpretedString ( obj )
%
%   Reads header interpretations from a comma separated values (.csv)
%   spreadsheet
%
%   Assumptions:
%       - Spreadsheet contains 2 columns
%       - First row of spreadsheet contains column names
%       - Column names are:
%         code, description
%

val='you fail';

%Which spreadsheet file has the definitions
definitionsFileName = [obj.headerType '.' obj.headerWord '.csv'];

%Does that spreadsheet file exist
if(exist(definitionsFileName,'file'))

    definitionsFile = File(definitionsFileName);

    %Get first row with column names, just in case (not needed?)
    k = textscan(definitionsFile.fid,'%q%q',1,'Delimiter',',');

    %Read the rest of the definitions
    v = textscan(definitionsFile.fid,'%q%q','Delimiter',',');
    M = containers.Map(v{1},v{2});

    try
        val = M(obj.value);
    catch ex
        %disp (ex.message)
        val = toString(obj);
    end    
else
    %return number
    val = toString(obj);    
end

return


%Get first row with column names, just in case (not needed if cols arre 
%in the right order)
%k = textscan(definitionsFile.fid,'%q%q',1,'Delimiter',',');

%Read the rest of the definitions
%v = textscan(definitionsFile.fid,'%q%q','Delimiter',',');


% [m n] = size(v)
% if m < 2
%   if strcmp(char(v(1)),code)
%       val = v(2);
%   else
%       val = code;
%   end
%   
%   return
% end

% rows = cat(2,v{:});
% 
% M = containers.Map(...
%     rows(1,:),...
%     rows(2,:));
% 
% code
% 
% val = codeRows(code);



end



