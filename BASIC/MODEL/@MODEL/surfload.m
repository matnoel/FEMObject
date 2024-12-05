function f = surfload(S,D,compo,fun,varargin)
% function f = surfload(S,D,compo,fun)
% assemblage d'une force surfacique
% D zone d'application de l'effort :
%       objet geometrique de meme dimension que la frontiere de S
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
fun = fcnchk(fun); % transformation de la fonction en inline (fonction matlab)

S0 = create_boundary(S);
if ~isempty(D)
    S0 = intersect(S0,D); % partie du maillage ou on applique l'effort
    S0 = createddlnode(S0,S.node); % on copie les ddl de S
end

f = calc_vector(S0,@load,compo,fun,varargin{:});

% try
P = calc_P(S,S0);
f = P'*f;
% catch
%     rep = findddl(S,'all',S0);
%     ftemp = f;
%     f = sparse(getnbddl(S),size(ftemp,2));
%     f(rep,:) = ftemp;
% end

if ~nofree
    f = freevector(S,f);
end
