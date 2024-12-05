function x = std(x)

Ex = expect(x);
x = expectnodimtimes([],x,x);
x = sqrt(simplify(x)-Ex.^2);

