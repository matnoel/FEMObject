function [ls,RVsample] = random(ls,varargin)
% function [ls,RVsample] = random(ls,varargin)

if ~israndom(ls)
    error('la LEVELSET n''est pas aleatoire')
else
    if iseval(ls)
        [ls.value,RVsample] = random(ls.value,varargin{:});
    else
        RV = RANDVARS(ls);
        
        RVsample = random(RV,varargin{:});
        ls = randomeval(ls,RVsample,RV);
    end
end
