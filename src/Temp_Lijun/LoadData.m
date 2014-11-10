%% Housing keeping
clear all
close all
clc

%% load data 1D
label = {'CCRB','EADB','FROB','GHIB','JCNB','JCSB','LCCB','MMNB','RMNB','SCYB','SMNB','VARB','VCAB'};


for i = 1 : length(label)
    [time,data,log] = fget_sac(['./data/BP.',label{i},'.BP1.SAC.bp']);
    T1D(i,:) = time;
    X1D(i,:) = data;
    info(i) = log;
end

% save data as MAT
save('SeismicData1D','T1D','X1D','info');

%% load data 1D
label = {'CCRB','EADB','FROB','GHIB','JCNB','JCSB','LCCB','MMNB','RMNB','SCYB','SMNB','VARB','VCAB'};


for i = 1 : length(label)
    for j = 1 : 3
        [time,data,log] = fget_sac(['./data/BP.',label{i},'.BP',num2str(j),'.SAC.bp']);
        index = (i-1)*3+j;
        T3D(index,:) = time;
        X3D(index,:) = data;
        info3D(index) = log;
    end
end

% save data as MAT
save('SeismicData3D','T3D','X3D','info3D');