function TT = sqSolver(a,b,c,rt) 
% solve the one variable square equation for Eikonal 3D
% a, b, c          see equation (2.3) in reference
% rt               righ-hand-side term: S(i,j)*dx

% Written by Entao Liu (liuentao@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology
%
%Reference: H. Zhao, A fast sweeping method for Eikonal equations,
%Mathematics of computation, 74(250), pp. 603-627, 2004

A=sort([a,b,c],'ascend');
if (A(2)-A(1)>=rt)
    TT=A(1)+rt;
elseif (0.5*(A(1)+A(2)+sqrt(2*rt^2-(A(1)-A(2))^2))<=A(3))
    TT=0.5*(A(1)+A(2)+sqrt(2*rt^2-(A(1)-A(2))^2));
else
    TT=(1/3)*(A(1)+A(2)+A(3)+sqrt((sum(A)^2-3*(A(1)^2+A(2)^2+A(3)^2-rt^2))));
end
