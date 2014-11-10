function obj = readDefinitions ( obj )
%
%function obj = readDefinitions (obj)
%
%   Reads header definitions from a comma separated values (.csv)
%   spreadsheet
%
%   Assumptions:
%       - Spreadsheet contains 8 columns
%       - First row of spreadsheet contains column names
%

k = textscan(obj.fid,'%q%q%q%q%q%q%q%q',1,'Delimiter',',');
obj.keys = cat(2,k{:});

v = textscan(obj.fid,'%q%q%q%q%q%q%q%q','Delimiter',',');
obj.values = cat(2,v{:});

end



