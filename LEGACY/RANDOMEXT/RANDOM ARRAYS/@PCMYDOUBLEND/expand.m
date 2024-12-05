function u=expand(u)

if isradial(u)
a = MULTIMATRIX(u.V,u.stodim);
s = size(a);
a = double(a)*double(u.L);

a = MULTIMATRIX(a,s,[size(a,2),1]);
u = PCMYDOUBLEND(a,getPC(u.L),u.stodim);

end
