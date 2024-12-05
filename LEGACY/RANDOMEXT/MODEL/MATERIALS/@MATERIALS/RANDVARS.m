function rvs = RANDVARS(mat)
% function rvs = RANDVARS(mat)

rvs = RANDVARS();
for k=1:mat.n
    if israndom(mat.MAT{k})
        rvs = RANDVARS(rvs,RANDVARS(mat.MAT{k}));
    end
end
