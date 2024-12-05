function ls=randomset(ls,varargin)
% function ls=randomset(ls,RVvalue,RV)
% remplace les valeurs des variables aleatoires pour les levelsetparametrees
% RVvalue : cell array (tableau de cellule des valeurs a affectees pour les variables)
% RV : indique les dimensions stochastiques de RVvalue
% RV peut etre un double de la meme taille que RVvalue ou encore un objet
% contenant les dimensions stochastiques (RANDVARS, POLYCHAOS, ...)

for k=1:ls.n
    ls.LS{k}=randomset(ls.LS{k},varargin{:});
end

