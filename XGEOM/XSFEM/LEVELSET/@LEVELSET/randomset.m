function ls=randomset(ls,varargin)
% function ls=randomset(ls,RVvalue,RV)
% remplace les valeurs des variables aleatoires pour les levelsetparametrees
% RVvalue : cell array (tableau de cellule des valeurs a affectees pour les variables)
% RV : indique les dimensions stochastiques de RVvalue
% RV peut etre un double de la meme taille que RVvalue ou encore un objet
% contenant les dimensions stochastiques (RANDVARS, POLYCHAOS, ...)

if ~isalevelset(ls)
 if isa(ls.value{1,1},'function_handle')
      levels = ls.value{1,2};
      for k=1:length(levels)
      levels{k} = randomset(levels{k},varargin{:}); 
      end
      ls.value{1,2}=levels ;
 else
      ls.value = randomsetparam(ls.value,varargin{:});
 end
 
end
    