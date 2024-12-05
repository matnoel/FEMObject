function rvs = RANDVARS(mat)
% function rvs = RANDVARS(mat)

rvs = RANDVARS();
for k=1:size(mat.param,1)
    if israndom(mat.param{k,2})
        rvs= RANDVARS(rvs,RANDVARS(mat.param{k,2}));
    end
end
