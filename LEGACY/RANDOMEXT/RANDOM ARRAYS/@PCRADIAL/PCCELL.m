function x = PCCELL(A)

V = A.V;
s=size(V{1});
L = A.D*double(vertcat(A.L{:}));
PC = getPC(A);
value = getmatrix(V);

value = value * sparse(L) ;
value = setmatrix(value,s(1),s(2));

x=PCCELL(value,PC);