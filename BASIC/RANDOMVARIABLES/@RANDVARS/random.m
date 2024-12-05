function A = random(u,varargin)
% function A = random(u,varargin)

n = getclassin('double',varargin,{1});
if ~isa(n,'cell')
    n = {n};
end

if ischarin('init',varargin);
    initstate
end

if ~isclassin('CONDRANDVAR',u.RV);
    A = cell(1,u.M);
    for k=1:u.M
        A{k} = random(u.RV{k},n{:});
    end
else % avec variables conditionnelles
    A = cell(1,u.M);
    [dummy,repcond] = isclassin('CONDRANDVAR',u.RV);
    [dummy,reprand] = isclassin('RANDVAR',u.RV);
    for k=1:length(reprand)
        A{reprand(k)} = random(u.RV{reprand(k)},n{:});
    end
    
    while ~isempty(repcond)
        for k=1:length(repcond)
            A{repcond(k)} = randomknowing(u.RV{repcond(k)},A(reprand),RANDVARS(u.RV(reprand)),n{:});
            if ~israndom(A{repcond(k)})
                reprand = union(reprand,repcond(k));
            end
        end
        repcond = setdiff(repcond,reprand);
    end
end
% A = [A{:}];
