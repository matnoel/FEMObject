function A = lhsrandom(u,n,p)

if nargin==2
    p=1;
end

A = lhsdesign(n,p);

A(:) = transfer(RVUNIFORM(0,1),u,A(:));

