function Z=expand(apcr)

apcr=doubleV(apcr);
V = apcr.V;
L = apcr.D*double(vertcat(apcr.L{:}));

n=size(V{1});
n = n(find(n>1));
if isempty(n)
    n=1;
end

l=size(L);

Z=zeros(prod(n),l(2));
for i=1:apcr.m 
  t = V{i};
  Z = Z + t(:)*L(i,:);
end

Z=reshape(Z,[n,l(2)]);

Z = PCARRAY(Z,apcr.POLYCHAOS);