function c = cellsum(a,k)
if ~isa(a,'cell')
    error('l''argument doit etre une cell')
end

sa = size(a);
if any(sa==1) & nargin==1
    a = a(:);k=1;
    sa = size(a);
elseif nargin==1
    k=1;
end

sc=sa;
sc(k)=1;
switch k
case 1
    c = cell(1,sa(2));     
    for j=1:sa(2)
        c{1,j} = a{1,j};
        for i=2:sa(1)
            c{1,j} = c{1,j} + a{i,j};   
        end
    end
case 2
    c = cell(sa(1),1); 
    for i=1:sa(1)
        c{i,1} = a{i,1};
        for j=2:sa(2)
            c{i,1} = c{i,1} + a{i,j};   
        end
    end
end



