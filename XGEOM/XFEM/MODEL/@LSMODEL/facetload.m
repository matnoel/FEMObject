function f = facetload(S,i,compo,fun,varargin)
% function f = surfload(S,i,compo,fun)
% assemblage d'une force surfacique
% i : numero de la facet
% compo : nom des composantes de la force
% fun : pointeur sur une fonction donnant la valeur des composantes en
%       colonnes (on peut avoir N colonnes correspondant a N cas de
%       chargement)
%       -> appel de fun(x,varargin{:}) ou x est un ND array
%       contenant les coordonnees de points ou on veut evaluer la fonction
%       size(x,1)=1
%       size(x,2)=dimension de l'espace de definition des points
%       size(x,3)*size(x,4)... : nombre de points
%       la sortie doit etre un ND array de taille length(compo)-by-N-by-size(x,3)-by-size(x,4)-by...

nofree = ischarin('nofree',varargin);

if isa(fun,'double')
    varargin{1} = fun;
    nargin = 1;
    % fun = inline('repmat(a,[1 1 size(x,3) size(x,4) size(x,5)])','x','a');
    fun = @(x,a) repmat(a,[1 1 size(x,3) size(x,4) size(x,5)]);
end
fun = fcnchk(fun); % tranformation de la fonction en inline (fonction matlab)

F = getfacet(S,i);
F = setmaterial(F,getmaterial(getgroupelem(S,1)));
F = actualise_ddl(F);
F = calcnumddl(F);
% construction de la structure des ddl des elements de F

if getnblevelsets(F)>0
    f = calc_vector(F,@loadls,getlevelsets(S),compo,fun,varargin{:});
else
    f = calc_vector(F,@load,compo,fun,varargin{:});
end

if ~nofree
    f = freevector(S,f);
end
