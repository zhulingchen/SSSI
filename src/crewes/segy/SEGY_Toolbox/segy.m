function [th bh xh]=segy(filename,varargin)
th=TraceHeader(filename);
bh=BinaryHeader(filename);
xh=TextHeader(filename);
%th=[];
end
