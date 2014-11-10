% Load the fbpick and offset file.
fprintf(1,'Loading fbpicks\n');
load fbpicks.txt;
fprintf(1,'Loading offset\n');
load offset.txt;
fprintf(1,'Done loading files\n');
% Make three vectors: trace#, offset, and fbpick
trace_vec = fbpicks(:,1)';
offset_vec = offset(:,2)';
fbpick_vec = fbpicks(:,2)';
% Geometry is 189 shots of 200 receivers
recs = 200;
shots = 189;
% Allocate a matrix of shots x receivers
shotoffset = zeros(shots, recs);
shotpick = zeros(shots, recs);
%display fbbreak in function of line location
load sinx.txt;
load siny.txt;
sin=sinx(:,1);
x=sinx(:,2)-sinx(189,2);
y=siny(:,2)-siny(189,2);
r=sqrt(x.^2+y.^2);
load uphole.txt;
u=uphole(:,2);
for i=1:shots
   first = (i-1)*recs+1;
   last = i*recs;
   fprintf(1,'Getting shot %d from traces %d to %d\n',i,first,last);
   shotoffset(i,1:recs) = offset_vec( first:last );
   shotpick(i,1:recs) = fbpick_vec( first:last );
   pickuphole(i,1:recs) = fbpick_vec( first:last) + u(i);
   reclocation(i,1:recs) = offset_vec( first:last ) + r(i);
end
