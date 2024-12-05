function rvs = RANDVARS(ls)
% function rvs = RANDVARS(ls)

rvs = RANDVARS();
if iseval(ls)
    if israndom(ls.value)
        rvs = RANDVARS(ls.value);
    end
else
    if isa(ls.value{1,1},'function_handle')
        levels = ls.value{1,2};
        for k=1:length(levels)
            rvs= RANDVARS(rvs,RANDVARS(levels{k}));
        end
    else
        
        for k=1:size(ls.value,1)
            if israndom(ls.value{k,2})
                rvs= RANDVARS(rvs,RANDVARS(ls.value{k,2}));
            end
        end
    end
end
