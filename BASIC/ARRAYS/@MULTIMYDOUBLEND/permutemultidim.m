function u = permutemultidim(u,k)
if k~=u.multidim
s=size(u.value);
p = 1:max([length(s),k,u.multidim]);
p(k)=u.multidim;
p(u.multidim)=k;
u.value = permute(u.value,p);
u.multidim=k;
end