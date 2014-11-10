function a=cell2char(c)
%
% a=cell2char(c)
% 
%Given a cell array of strings, cell2char converts these into a single
%string vector with each cell string separated by a logical '|'. For example,
%let c{1}='fred' and c{2}='billygoat'. Then 
% a=cell2char(c) results in a='fred|billygoat'. The final form is usedful
% in user interface dialogs like askthings.
%

n=length(c);

nmax=n*120;%assume 120 charaters in each cell

a=char(32*ones(1,nmax));
n1=1;
for k=1:n
    sk=c{k};%the kth string
    if(k<n)
        n2=n1+length(sk);
        a(n1:n2)=[sk '|'];
    else
        n2=n1+length(sk)-1;
        a(n1:n2)=sk;
    end
    n1=n2+1;
end

if(n1<length(a))
    a(n1:end)='';
end