
[out_x, out_y] = gradient(buffer);
out_x = out_x(:);
out_y = out_y(:);
A = [out_x,out_y]';

AA = A*A';

a = AA(1);
b = AA(2);
c = AA(4);

d = sqrt((a+c)^2 - 4*(a*c - (b)^2));

lam_1 = (a+c-d)/2;
lam_2 = (a+c+d)/2;

eig_1(1) = -b;
eig_1(2) = (a-lam_1); eig_1 = eig_1/sqrt(eig_1*eig_1');

eig_2(1) = -b;
eig_2(2) = (a-lam_2); eig_2 = eig_2/sqrt(eig_2*eig_2');

theta = -atan(eig_1(2)/eig_1(1));

R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

H = R*S*(S')*(R');


for l_x = [-50:50]
 for l_y = [-50:50]
    
    l
x = x_loc + l_x;
y = y_loc + l_y;

H = [C00(y_loc+y_shift, x_loc+x_shift), C01(y_loc+y_shift, x_loc+x_shift); C01(y_loc+y_shift, x_loc+x_shift), C11(y_loc+y_shift, x_loc+x_shift)];
X = [x - (x_loc + x_shift), y - (y_loc + y_shift)]';
h = 4;
weight(l_y+51, l_x+51) = exp(-1/(2*h^2)*X'*inv(H)*X);
 end
end

