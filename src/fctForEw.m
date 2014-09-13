function D_R = fctForEw(ID0,ID1)
% This function is implementing FCT (flux-correct transport) to eliminate
% the unmerical dispersion for Elastic wavefield.
% This function is just for Two-Dimension or Three-Dimension 
% Strigged-Grid Difference
% Author: Yan Xinfei(fulcrum)
% Finishing Time:2007/10/20
% Usage:
%       D_R = FCTforEW(ID0,ID1) ID0 and ID1 are the data of (n-1/2)th 
%            and (n+1/2)th layer respectively.
d = ndims(ID0);% The Dimension of Matrix
[m,n] = size(ID0);% The Size of Matrix
if nargin ~= 2
    error('The number of input arguments is wrong!');
else
    if d == 2
        % Step 1 & Step 2
        et1 = 0.03;et2 = 0.00;% values between 0.02 and 0.1
        Qx0 = zeros(m-1,n);Qy0 = zeros(m,n-1);
        Qx1 = zeros(m-1,n);Qy1 = zeros(m,n-1);
        for i = 1:m
            for j= 1:n
                if i~=m && j~=n
                   Qx0(i,j) = et1*(ID0(i+1,j)-ID0(i,j));
                   Qx1(i,j) = et2*(ID1(i+1,j)-ID1(i,j));
                   Qy0(i,j) = et1*(ID0(i,j+1)-ID0(i,j));
                   Qy1(i,j) = et2*(ID1(i,j+1)-ID1(i,j));
               elseif i==m && j~=n
                   Qy0(i,j) = et1*(ID0(i,j+1)-ID0(i,j));
                   Qy1(i,j) = et2*(ID1(i,j+1)-ID1(i,j));
               elseif j==n && i~=m
                   Qx0(i,j) = et1*(ID0(i+1,j)-ID0(i,j));
                   Qx1(i,j) = et2*(ID1(i+1,j)-ID1(i,j));
               end
            end
        end
        
        % Step 3
        ID_1 = zeros(m,n);
        for i = 2:m-1
            for j = 2:n-1
                ID_1(i,j) = ID1(i,j) + (Qx0(i,j)-Qx0(i-1,j)) + (Qy0(i,j)-Qy0(i,j-1));
            end
        end
        
        % Step 4
        X = zeros(m-1,n);Y = zeros(m,n-1);
        for i = 1:m
            for j = 1:n
                if i~=m && j~=n
                   X(i,j) = ID_1(i+1,j)-ID_1(i,j);
                   Y(i,j) = ID_1(i,j+1)-ID_1(i,j);
               elseif i==m && j~=n
                   Y(i,j) = ID_1(i,j+1)-ID_1(i,j);
               elseif j==n && i~=m
                   X(i,j) = ID_1(i+1,j)-ID_1(i,j);
               end
            end
        end
        
        %Step 5
        D_R = ID1;
        for i = 3:m-2
            for j = 3:n-2
                Sx0 = sign(Qx1(i-1,j));Sx1 = sign(Qx1(i,j));
                Sy0 = sign(Qy1(i,j-1));Sy1 = sign(Qy1(i,j));
                Xc0=Sx0*max([0,min([Sx0*X(i-2,j),abs(Qx1(i-1,j)),Sx0*X(i,j)])]);
                Xc1=Sx1*max([0,min([Sx1*X(i-1,j),abs(Qx1(i,j)),Sx1*X(i+1,j)])]);
                Yc0=Sy0*max([0,min([Sy0*Y(i,j-2),abs(Qy1(i,j-1)),Sy0*Y(i,j)])]);
                Yc1=Sy1*max([0,min([Sy1*Y(i,j-1),abs(Qy1(i,j)),Sy1*Y(i,j+1)])]);
                D_R(i,j)=ID_1(i,j)-(Xc1-Xc0)-(Yc1-Yc0);
            end
        end
    end
end