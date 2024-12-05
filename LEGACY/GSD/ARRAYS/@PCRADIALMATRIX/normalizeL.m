function u = normalizeL(u)

D = sparse(u.m,u.m);
for i=1:u.m
D(i,i) = norm(u.L(i));
end

u.L = inv(D) * u.L ; 
u.V = MULTIMATRIX(double(u.V)*D,size(u),sizem(u.V));
u.DLmasse={};

 