function M = reduce(A,j)

Z = A.F{j};
%A.F(j) = [];
a = A.alpha;
n = size(a);
for i=1:A.dim
    if i~=j
a = ttv(a,full([A.F{i}{:}]'),i);
n(i)=1;
a = reshape(a,n);
    end
end
a = double(squeeze(a));

global mexcompiled;

if 1%~mexcompiled
    M = a(1)*Z{1};
    for i=2:size(a,1)
        M = M + a(i)*Z{i};
    end
else
    M = redux(a,Z);
end


return
