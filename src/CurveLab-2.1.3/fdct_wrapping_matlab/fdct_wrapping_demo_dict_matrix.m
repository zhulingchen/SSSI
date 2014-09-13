close all;
clear;
clc;

load velocityModelMarmousiTest.mat
[nx, ny] = size(velocityModel);
is_real = 1;

%forward curvelet transform
disp('Take curvelet transform: fdct_wrapping');
tic;
C = fdct_wrapping(velocityModel, is_real);
toc;

% vectorize curvelet coefficients
nc = 0;
vc = [];
for s=1:length(C)
    for w=1:length(C{s})
        nc = nc + length(abs(C{s}{w}(:)));
        vc = [vc; C{s}{w}(:)];
    end
end

C0 = C;
for s=1:length(C0)
    for w=1:length(C0{s})
        [A, B] = size(C0{s}{w});
        C0{s}{w} = zeros(A, B);
    end
end

Phi = zeros(nx * ny, nc);
icol = 1;
for s=1:length(C0)
    for w=1:length(C0{s})
        [A, B] = size(C0{s}{w});
        % pay attention: rows first to fit for the order used in C{s}{w}(:)
        for b = 1:B
            for a = 1:A
                display(icol);
                C0{s}{w}(a, b) = 1;
                tmp = ifdct_wrapping(C0, is_real, nx, ny);
                Phi(:, icol) = tmp(:);
                icol = icol + 1;
                C0{s}{w} = zeros(A, B);
            end
        end
    end
end

vm2 = Phi * vc;
vm2 = reshape(vm2, nx, ny);
delta = velocityModel - vm2;
display(max(abs(delta(:))));

img = fdct_wrapping_dispcoef(C);
figure; imagesc(abs(img));


disp('Take inverse curvelet transform: ifdct_wrapping');
tic;
velocityModelRec = ifdct_wrapping(C, is_real, nx, ny);
toc;

delta = velocityModel - velocityModelRec;
max(abs(delta(:)))

