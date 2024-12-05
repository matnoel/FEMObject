function ls = randomeval(ls,varargin)
% function ls = randomeval(ls,x,RV)
% calcul de realisations des LEVELSETS
% on applique randomeval(.,x,RV) a toutes les LEVELSET
% 
%  See also LEVELSET/randomeval

for k=1:ls.n
    ls.LS{k} = randomeval(ls.LS{k},varargin{:});
end

