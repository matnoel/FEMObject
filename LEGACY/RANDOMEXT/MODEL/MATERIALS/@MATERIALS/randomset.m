function mat = randomset(mat,varargin)
% function mat = randomset(mat,RVvalue,RV)
% remplace les valeurs des RANDVAR pour les material parametrees
% RVvalue : cell array (tableau de cellule des valeurs a affectees pour les variables)
% RV : indique les dimensions stochastiques de RVvalue
% RV peut etre un double de la meme taille que RVvalue ou encore un objet
% contenant les dimensions stochastiques (RANDVARS, POLYCHAOS, ...)
%
% See also MATERIAL/randomset

for k=1:mat.n
    mat.MAT{k} = randomset(mat.MAT{k},varargin{:});
end

