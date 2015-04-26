% SETPATH adds necessary directories to search path
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

homeDir = fileparts(pwd);
addpath(genpath([homeDir, '/model_data']));
addpath(genpath([homeDir, '/src']));
if ~isunix
    rmpath(genpath([homeDir, '/src/CurveLab-2.1.3/fdct_usfft_cpp']));
    rmpath(genpath([homeDir, '/src/CurveLab-2.1.3/fdct_wrapping_cpp']));
    rmpath(genpath([homeDir, '/src/CurveLab-2.1.3/fdct3d']));
end

model_data_path = [homeDir, '/model_data'];
images_path = [homeDir, '/images'];
videos_path = [homeDir, '/videos'];
