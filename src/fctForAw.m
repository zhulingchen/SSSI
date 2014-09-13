function D_R = fctForAw(data0,data1,data2)
% FCTFORAW is implementing FCT (flux-correct transport) to eliminate
% the unmerical dispersion for Acoustic wavefield. 
% This function is just for Two-Dimension or Three-Dimension 
% Regular-Grid Difference
%
% Author: Yan Xinfei(fulcrum)
% Finishing Time: 2007/10/20
%
% Usage:
%       D_R = FCTforAW(data0,data1,data2) data0 is the data of (n-1)th
%       layer,and data1 and data2 are the the data of (n)th and (n+1)th layer respectively.

d = ndims(data2);% The Dimension of Matrix
[m,n] = size(data2);% The Size of Matrix
if nargin ~= 3
    error('The number of input arguments is wrong!');
else
    if d == 2
        % Step 1 & Step 2
        et1=0.04;et2=0.045;% values between 0.02 and 0.1
        Qxn=zeros(m-1,n);Qyn=zeros(m,n-1);Q_xn=Qxn;Q_yn=Qyn;Xt=Qxn;Yt=Qyn;
        data_t=zeros(m,n);data_r=data_t;
        for i = 1:m
            for j= 1:n
                if i~=m && j~=n
                   Qxn(i,j)=et1*((data1(i+1,j)-data1(i,j))-(data0(i+1,j)-data0(i,j)));
                   Qyn(i,j)=et1*((data1(i,j+1)-data1(i,j))-(data0(i,j+1)-data0(i,j)));
                   Q_xn(i,j)=et2*((data2(i+1,j)-data2(i,j))-(data1(i+1,j)-data1(i,j)));
                   Q_yn(i,j)=et2*((data2(i,j+1)-data2(i,j))-(data1(i,j+1)-data1(i,j)));
               elseif i==m && j~=n
                   Qyn(i,j)=et1*((data1(i,j+1)-data1(i,j))-(data0(i,j+1)-data0(i,j)));
                   Q_yn(i,j)=et2*((data2(i,j+1)-data2(i,j))-(data1(i,j+1)-data1(i,j)));
               elseif j==n && i~=m
                   Qxn(i,j)=et1*((data1(i+1,j)-data1(i,j))-(data0(i+1,j)-data0(i,j)));
                   Q_xn(i,j)=et2*((data2(i+1,j)-data2(i,j))-(data1(i+1,j)-data1(i,j)));
               end
            end
        end
        
        % Step 3
        for i = 2:m-1
            for j = 2:n-1
                data_t(i,j)=data2(i,j)+(Qxn(i,j)-Qxn(i-1,j))+(Qyn(i,j)-Qyn(i,j-1));
            end
        end
        
        % Step 4
        for i = 1:m
            for j = 1:n
                if i~=m && j~=n
                   Xt(i,j)=(data_t(i+1,j)-data_t(i,j))-(data_r(i+1,j)-data_r(i,j));
                   Yt(i,j)=(data_t(i,j+1)-data_t(i,j))-(data_r(i,j+1)-data_r(i,j));
               elseif i==m && j~=n
                   Yt(i,j)=(data_t(i,j+1)-data_t(i,j))-(data_r(i,j+1)-data_r(i,j));
               elseif j==n && i~=m
                   Xt(i,j)=(data_t(i+1,j)-data_t(i,j))-(data_r(i+1,j)-data_r(i,j));
               end
            end
        end
        
        %Step 5
        D_R=data2;
        for i = 3:m-2
            for j = 3:n-2
                Sx0=sign(Q_xn(i-1,j));Sx1=sign(Q_xn(i,j));Sy0=sign(Q_yn(i,j-1));Sy1=sign(Q_yn(i,j));
                Xc0=Sx0*max([0,min([Sx0*Xt(i-2,j),abs(Q_xn(i-1,j)),Sx0*Xt(i,j)])]);
                Xc1=Sx1*max([0,min([Sx1*Xt(i-1,j),abs(Q_xn(i,j)),Sx1*Xt(i+1,j)])]);
                Yc0=Sy0*max([0,min([Sy0*Yt(i,j-2),abs(Q_yn(i,j-1)),Sy0*Yt(i,j)])]);
                Yc1=Sy1*max([0,min([Sy1*Yt(i,j-1),abs(Q_yn(i,j)),Sy1*Yt(i,j+1)])]);
                D_R(i,j)=data_t(i,j)-(Xc1-Xc0)-(Yc1-Yc0);
            end
        end
        data_r=data_t;
    end
end