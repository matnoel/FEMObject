function A = lhsrandom(u,n,j)

A = lhsdesign(n,u.M);
A = transfer(RANDVARS(RVUNIFORM(0,1),u.M),u,A);
A = mat2cell(A,n,ones(1,u.M));

