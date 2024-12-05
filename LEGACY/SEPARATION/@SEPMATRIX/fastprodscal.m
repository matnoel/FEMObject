function ps = fastprodscal(u,v)
u=u';
ps=sum(prod([cell2mat(multiplyF(u.F,v.F)) kron(u.alpha(1:u.m),v.alpha(1:v.m))'],2));