function header = reshapeHeader( header, dimension )
%function obj = reshapeHeader( header, dimension )
%   Reshape text header
%  Where:
%   dimension = '1D' or '2D'
%

if(~ischar(header))
    error ('Textual header must be of type char');
end

[m n]=size(header);

switch(m+n)
    case(3201)     % 1D vector input
        if (n > m) % input is row vector
            header = header'; % always start with a column vector
        end
        switch(dimension)
            case('1D')
                   % all's well, do nothing else
            case('2D')
                header = reshape(header,80,40)';
            otherwise
                error('dimension must be ''1D'' or ''2D''');
        end
    case(120)      % input is 2D matrix
        if (m > n) % input is rotated?
            error('Input must be 40 rows by 80 columns');
        end
        switch(dimension)
            case('1D')
                header = reshape(header',3200,1);
            case('2D')
                   % all's well, do nothing
            otherwise
                error('dimension must be either ''vector'' or ''matrix''');
        end
    otherwise
       error ('Textual header must be 3200 bytes');
end

end



