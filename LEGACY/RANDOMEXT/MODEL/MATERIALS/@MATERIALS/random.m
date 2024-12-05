function mat = random(mat,varargin)
% function mat = random(mat,varargin)

if ~israndom(mat)
    error('la MATERIALS n''est pas aleatoire')
else
    RV = RANDVARS(mat);
    RVsample = random(RV,varargin{:});
    for k=1:mat.n
        mat.MAT{k} = randomeval(mat.MAT{k},RVsample,RV);
    end
end
