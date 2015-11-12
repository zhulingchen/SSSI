function dirClass = getPatchGradClass(yin, dirRange)
% GETPATCHGRADCLASS Estimate gradient class of a patch
% GETPATCHGRADCLASS estimates the class of image gradient of the input
% patch
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

[gx, gy] = gradient(yin);
gDir = atand(gx ./ gy);
binCount = histc(gDir(:), dirRange);
binCount = binCount(1:end-1);
binMiddle = length(binCount)/2 + 1;
binCount_new = binCount(1:binMiddle - 1) + binCount(end:-1:binMiddle);
[~, dirClass] = max(binCount_new);