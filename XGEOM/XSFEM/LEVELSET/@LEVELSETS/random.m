function ls = random(ls,varargin)
% function ls = random(ls,varargin)

if ~israndom(ls)
    error('la LEVELSETS n''est pas aleatoire')
else
    RV = RANDVARS(ls);
    x = random(RV,varargin{:});
    for k=1:ls.n
        ls.LS{k}=randomeval(ls.LS{k},x,RV);
    end
end
