function shot = shotrange(string)
[s1, rest] = strtok(string,',');
[s2, rest] = strtok(rest,',');
a = sscanf(s1, '(%d:%d)');
if(isempty(a))
   a = sscanf(s1, '%d');
   i1 = a;
   i2 = a;
else
   i1 = a(1);
   i2 = a(2);
end
b = sscanf(s2, '(%d:%d)');
if(isempty(b))
   b = sscanf(s2, '%d');
   j1 = b;
   j2 = b;
else
   j1 = b(1);
   j2 = b(2);
end
% Now, 'a' contains the range of the first shots (either [1] or [1 2]),
% and 'b' contains the range of the 2nd shots.
% Generate a list of shot pairs
k=0;
for i=i1:i2
   for j=j1:j2
      if( i<j )
         k = k+1;
         shot(k,1) = i;
         shot(k,2) = j;
      end
   end
end
