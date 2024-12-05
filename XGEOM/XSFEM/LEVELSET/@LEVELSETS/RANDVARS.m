function rvs = RANDVARS(ls)
% function rvs = RANDVARS(ls)

rvs = RANDVARS();
for k=1:ls.n
    if israndom(ls.LS{k})
        rvs = RANDVARS(rvs,RANDVARS(ls.LS{k}));
    end
end
