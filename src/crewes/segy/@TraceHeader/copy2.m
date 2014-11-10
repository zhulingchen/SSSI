function objout=copy2(objin)
% this function will make a value copy of the traceheader class;

objout=TraceHeader(objin.filename,'new','1');
objout.traceoffsets=objin.traceoffsets;
objout.tracetype=objin.tracetype;
objout.header=objin.header;
objout.definitions=objin.definitions;
objout.hdrsize=objin.hdrsize;
objout.hdroffset=objin.hdroffset;
objout.nontypecasthdr=objin.nontypecasthdr;
objout.filefmt=objin.filefmt;
objout.machineformat=objin.machineformat;
end