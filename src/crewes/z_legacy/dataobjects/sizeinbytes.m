function bytes=sizeinbytes(thing)

klass=class(thing);

bytes=8*prod(size(thing));%default size

if(strcmp(klass,'char'))
   bytes=2*prod(size(thing));
elseif(strcmp(klass,'int8')|strcmp(klass,'uint8'))
   bytes=prod(size('thing'));
elseif(strcmp(klass,'uint16')|strcmp(klass,'int16'))
   bytes=2*prod(size(thing));
elseif(strcmp(klass,'int32')|strcmp(klass,'uint32'))
   bytes=4*prod(size(thing));
elseif(strcmp(klass,'struct'))
    [l,m]=size(thing);
    names=fieldnames(thing);
    numnames=size(names,1);
    bytes=0;
    for k=1:numnames
        bytes=bytes+sizeinbytes(getfield(thing,char(names(k,:))));
    end
 end
