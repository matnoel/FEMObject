function Z=getpccompo(apcr,p)

apcr=doubleV(apcr);
V = apcr.V;
L = getpccompo(vertcat(apcr.L{:}),p);
L = apcr.D*L;

n=size(V{1});
l=size(L);

Z=sparse(prod(n),l(2));
for i=1:apcr.m 
  t = V{i};
  Z = Z + t(:)*sparse(L(i,:));
end

Z=reshape(Z,[n,l(2)]);
Z=squeeze(Z);
