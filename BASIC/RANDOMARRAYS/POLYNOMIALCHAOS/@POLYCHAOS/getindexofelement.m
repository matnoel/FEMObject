function ind = getindexofelement(PC,e)

H = PC.RANDPOLYS;
if PC.M>1
    error('pas programme')
end

h = getpoly(H,1);
indices = PC.indices(:,1);
n = getnbelem(h);
ind = find(mod(indices,n)+1==e);



