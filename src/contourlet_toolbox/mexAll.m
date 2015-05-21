close all;
clear;
clc;

cFilesStr = dir('*.c');
cFilenames = {cFilesStr.name};

for ii = 1:length(cFilenames)
    mex(cFilenames{ii});
end