function f = surfload_elem(S,D,compo,fun,PC,varargin)
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
%       la sortie doit etre un ND array de taille
%       length(compo)-by-N-by-size(x,3)-by-size(x,4)-by...

if isa(fun,'double')
    % varargin{1} = fun;
    % nargin = 1;
    varargin = [fun,varargin];
    % fun = inline('repmat(a,[1 1 size(x,3) size(x,4) size(x,5)])','x','a');
    fun = @(x,a) repmat(a,[1 1 size(x,3) size(x,4) size(x,5)]);
end
fun = fcnchk(fun); % tranformation de la fonction en inline (fonction matlab)

S0 = create_boundary(S);
S0 = intersect(S0,D); % partie du maillage ou on applique l'effort
% S0 = createddlnode(S0,S.node); % on copie les ddl de S
% keyboard
% tol = getcharin('tol',varargin,1e-1);
% varargin = delcharin('tol',varargin);
% S0 = final_elem(S0,'tolsplit',tol);
% keyboard
% S0 = final(S0); % à voir

if ~israndom(S0.ls)
    if getnblevelsets(S0)>0
        warning('attention l''enrichissenent LEVELSET n''est pas conservé sur le bord')
        f = calc_multivector(S0,@lsload,S0.ls,compo,fun,varargin{:});
    else
        f = calc_vector(S0,@load,compo,fun,varargin{:});
    end
else
    p = S0.nbgroupelem;
    if p==1
        switch getlstype(S0.groupelem{1})
            case {'in'}
                f = calc_vector(S0,@load,compo,fun,varargin{:});
            case {'cut'}
                f = calc_pcvector(S0,PC,@lsloadpc,S0,compo,fun,varargin{:});
            case {'out'}
                f = calc_vector(S0,@load,compo,fun,varargin{:});
        end
    else
        PC = getPC(PC);
        f = calc_pcvector(S0,PC,@lsloadpc,S0,compo,fun,varargin{:});
        % f = calc_pcvector_elem(S0,PC,@lsloadpc_elem,S0,compo,fun,varargin{:});
    end
end
% calc_pcvector à faire si S0.ls random
% f = calc_pcvector(S0,PC,@lsloadpc,S0.ls,compo,fun,varargin{:});
P = calc_P(S,S0);
f = P'*f;

f = freevector(S.BCOND,f);
