function fe = loadls(elem,node,ls,varargin)
% function fe = loadls(elem,node,ls,compo,fun,varargin)
% compo : nom des composantes de la force
% fun : pointeur sur une fonction donnant la valeur des composantes en
% colonnes (on peut avoir N colonnes correspondant a N cas de
%       chargement)
% appel de fun(x,varargin{:}) ou x sont le tableau de coordonnees de P points ou on veut evaluer la fonction
% size(x,1)=P size(x,2)=dimension de l'espace de definition des points
% la sortie est un 3D array de taille length(compo)*N*P  

if isempty(getlsnumber(elem))
    nature = 'domain';
    ls = [];  
else
    ls = getlevelset(ls,getlsnumber(elem));
    nature = getnature(ls);
end

switch nature
    case 'crack'
        if isenrich(elem)
            error('pas prevu')
        else
        fe = loadlscrack(elem,node,ls,varargin{:});
        end      
    case 'material'
        fe = loadlsmaterial(elem,node,ls,varargin{:});
    case 'domain'
        fe = loadlsdomain(elem,node,ls,varargin{:});
end

