function M = reduce(A,j)

if A.dim==1
    M = A.F{1}*A.alpha(1);
    for i=2:A.m
        M = M + A.F{i}*A.alpha(i);
    end
else
    Z = A.F(:,j);
    A.F(:,j) = [];
    a = A.alpha.*prod(cell2mat(A.F),2)';
    if 1%~mexcompiled
         M = a(1)*Z{1};
         for i=2:A.m
             M = M + a(i)*Z{i};
         end
    else
        M = redux(a,Z);
    end
end
return
