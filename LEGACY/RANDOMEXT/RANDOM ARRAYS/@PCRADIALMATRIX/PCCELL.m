function x = PCCELL(A)

V = A.V;
s=size(V);
L = A.D*double(A.L);
PC = getPC(A);

value = double(V) * sparse(L) ;
for k=1:length(PC)
x{k} = reshape(value,s);
end

x=PCCELL(x,PC);