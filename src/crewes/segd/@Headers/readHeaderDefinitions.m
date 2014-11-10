function obj = readHeaderDefinitions ( obj )
%
%function obj = readDefinitions (obj)
%
%   Reads header definitions from a comma separated values (.csv)
%   spreadsheet
%
%   Assumptions:
%       - Spreadsheet contains 7 columns
%       - First row of spreadsheet contains column names
%       - First column of spreadsheet contains header word names
%       - column names are:
%         name,startByte,endByte,startNibble,endNibble,format,description
%

%Which spreadsheet file has the definitions
definitionsFile = File([obj.headerType '.csv']);

%Get the column names (first row of spreadsheet)
k = textscan(definitionsFile.fid,'%q%q%q%q%q%q%q',1,'Delimiter',',');
obj.headerDefCols = containers.Map(cat(2,k{:}),1:length(k));

%Read the rest of the spreadsheet
v = textscan(definitionsFile.fid,'%q%q%q%q%q%q%q','Delimiter',',');
obj.headerDefs = cat(2,v{:});

%Get the row names (first column of spreadsheet)

[m n] = size(obj.headerDefs);

obj.headerDefRows = containers.Map(...
    obj.headerDefs(:,obj.headerDefCols('name')),...
    1:m);

definitionsFile.closeFile();

end



