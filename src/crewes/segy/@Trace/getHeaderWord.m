function word = getHeaderWord ( obj, wordname)
%word = getHeaderWord
%   Detailed explanation goes here

%get row index in definitions spreadsheet for wordname
w = strmatch(upper(wordname),upper(obj.definitions.values(:,obj.definitions.n)));

%get info for header word
sbyte = str2double(obj.definitions.values{w,obj.definitions.s});
ebyte = str2double(obj.definitions.values{w,obj.definitions.e});
type  = obj.definitions.values{w,obj.definitions.t};

%get header words
word  = obj.header(sbyte:ebyte,:);
[m n] = size(word);
word  = typecast( reshape(word,m*n,1),type );

%byte order
[c m e] = computer;

%assume segy is big-endian
if (strcmp(e,'L'))
    word=swapbytes(word);
end

end

