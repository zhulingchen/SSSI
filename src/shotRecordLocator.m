function [fileInd, colInd] = shotRecordLocator(shotNumber,shotArray,shotsInFile)

% Number of traces per shot record
nTraces = shotArray{1}.NumberOfTraces/shotsInFile(1);
ds = diff([0 shotsInFile]);

% Determine file record to return
fileInd = sum(shotsInFile < shotNumber);
if fileInd < length(ds)
    fileInd = fileInd + 1;
end

% Determine column indicies (traces) for shot record (stored in reverse
% order)

skip = ds(fileInd) - (shotsInFile(fileInd) - shotNumber);
colInd = (1:nTraces) + nTraces*(skip-1);




