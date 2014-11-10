function val=checkforrightbyteorder(val,filefmt)
if nargin<2
    filefmt='B';
end
%byte order
    [ab, ac, e] = computer;
 %assume segy is big-endian
        if (~strcmp(e,filefmt))
            val=swapbytes(val);
        end
end