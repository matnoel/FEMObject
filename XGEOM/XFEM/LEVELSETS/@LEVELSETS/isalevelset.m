function rep = isalevelset(ls)
warning('obsolete : utiliser iseval')
rep = 1; 

for k=1:ls.n
rep=rep & isalevelset(ls.LS{k});
end

