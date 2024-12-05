function ls = randomeval(ls,varargin)
% function ls = randomeval(ls,x,RV)
% calcul de realisations de la LEVELSET
% RV : objet  (RANDVARS, POLYCHAOS, ...) indiquant les dimensions stochastiques de x
% par defaut RV = RANDVARS(ls)
% x (n-by-M double ou 1-by-M cell de n-by-1 double) contenant les realisations des M variables RV
% - si ls est une levelset et ls.value est une PCMATRIX
%   on appelle randomeval(ls.value,x,RV) (realisations de la PCMATRIX)
% - sinon on applique randomeval(.,x,RV) a tous les parametres de la
%   levelset
%
%  See also LEVELSETS/randomeval, PCMATRIX/randomeval, RANDVAR/randomeval

if israndom(ls)
    if iseval(ls)
        ls.value = randomeval(ls.value,varargin{:});
    else
        if isa(ls.value{1,1},'function_handle')
            levels = ls.value{1,2};
            for k=1:length(levels)
                levels{k} = randomeval(levels{k},varargin{:});
            end
            ls.value{1,2} = levels ;
        else
            ls.value = funrandomparam(ls.value,@randomeval,varargin{:});
        end
        
    end
end
