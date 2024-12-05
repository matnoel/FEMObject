function nbelemsto = getnbelemsto(M)

if israndom(M)
nbelemsto = zeros(M.nbgroupelem,1);

for p=1:M.nbgroupelem
    temp = zeros(getnbelem(M.groupelem{p}),1);
   for j=1:getnbelem(M.groupelem{p})
       stomesh = getparam(getelem(M.groupelem{p},j),'stomesh');
       temp(j) = stomesh{1}.n;
   end
   nbelemsto(p) = sum(temp);
end
nbelemsto = sum(nbelemsto);
else
    error('le modèle n''est pas aléatoire')
end