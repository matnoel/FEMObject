function Z=expand(A)

s=size(A);
PC = getPC(A);
A = getmatrix(A.value);
if all(s>1);
A=reshape(full(A),[s size(A,2)]);
end

Z = PCARRAY(A,PC);
