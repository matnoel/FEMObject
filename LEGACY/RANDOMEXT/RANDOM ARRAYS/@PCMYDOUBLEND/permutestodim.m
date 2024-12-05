function u = permutestodim(u,k)
if k~=u.stodim
s=size(u.V);
p = 1:max([length(s),k,u.stodim]);
p(k)=u.stodim;
p(u.stodim)=k;
u.V = permute(u.V,p);
u.stodim=k;
end