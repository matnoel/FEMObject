function Z = integrate(a)

if isa(a.value,'cell')
error('pas programme')
else
M = getMmatrix(a);    
Z = reshape(sum(a.value*M,2),a.s);
end