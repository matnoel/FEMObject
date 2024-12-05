function f = bodyload(S,D,compo,fun,varargin)
% function f = bodyload(S,D,compo,fun,varargin)
% assemblage d'une force volumique
% D zone d'application de l'effort :
%       objet geometrique de meme dimension que S
% compo : nom des composantes de la force
% fun : pointeur sur une fonction donnant la valeur des composantes en
%       colonnes, peut N colonnes correspondant N cas de chargement
%       -> appel de fun(x,varargin{:}) ou x sont le tableau de coordonnees de P points
%       ou on veut evaluer la fonction
%       size(x,1)=P size(x,2)=dimension de l'espace de definition des points
%       la sortie est un 3D array de taille length(compo)*N*P

nofree = ischarin('nofree',varargin);

if isa(fun,'double')
    varargin{1} = fun;
    nargin = 1;
    % fun = inline('repmat(a,[1 1 size(x,3) size(x,4) size(x,5)])','x','a');
    fun = @(x,a) repmat(a,[1 1 size(x,3) size(x,4) size(x,5)]);
end
fun = fcnchk(fun); % transformation de la fonction en inline (fonction matlab)

if ~isempty(D)
    S0 = intersect(S,D); % partie du maillage ou on applique l'effort
    S0 = createddlnode(S0,S.node); % on copie les ddl de S
else
    S0 = S;
end
% nume est le numero des elements de S

fprintf('-> Compute element forces ... ')
for p=1:S0.nbgroupelem
    if length(S.ls)>0
        fe{p} = lsload(S0.groupelem{p},S.node,S,compo,fun,varargin{:});
    else
        fe{p} = load(S0.groupelem{p},S.node,S,compo,fun,varargin{:});
    end
    pourcentage(p,S0.nbgroupelem);
end

f = assemble_vectorelem(S,fe);

% f = FEVECTOR(f,S);
if ~nofree
    f = freevector(S.BCOND,f);
end
