load SeismicData3D.mat
x = squeeze(X3D(1,1:1000,:));

WL = 10; % samples
overlap = 0.5; % overlap percentage

p = [];
step = floor(WL*(1-overlap));
Nb = floor(length(x)/step-1);
s = 1:step:length(x)-step;

for i = 0 : Nb-1
    frame = x(1+i*step:WL+i*step,:);
    p = [p getpol(frame)];
end
figure(gcf);
subplot(2,1,1);
plot(x(:,1));
subplot(2,1,2);
stem(s,p);axis([0,250,0.5,1])