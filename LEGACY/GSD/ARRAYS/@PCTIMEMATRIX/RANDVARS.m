function RV = RANDVARS(u)
% function RV = RANDVARS(u)

if ~israndom(u)
    error('la timematrix n''est pas aleatoire')
else
    if isa(u.value,'cell')
        for i=1:length(u.value)
            if isa(u.value{i},'POLYCHAOS')
                RV= RANDVARS(u.value{i});
                break
            end
        end
    else
        RV  = RANDVARS(u.value);
    end
end
