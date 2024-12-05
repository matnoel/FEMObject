function A = PCMATRIX(R)

s=size(R.V{1});
for k=1:R.m
   V{k}=R.V{k}(:); 
end
V = [V{:}];
PC = R.POLYCHAOS ;
L = double(vertcat(R.L{:}));

A = V * R.D * L;

A = PCMATRIX(A,s,PC);