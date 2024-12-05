function h = POLYFE(n,p,type,varargin)
% function h = POLYFE(n,p)
% base element fini sur [0,1]
% polynomes orthonorm�s sur [0,1] pour la mesure 1
% p : ordre du polynorme (ne sert pas a grand chose)
% n : si scalaire : decoupage uniforme en n elements
%     si vecteur  : coordonnes des noeuds du maillage
%
% function h = POLYFE(n,p,'node')
% n correspond aux coordonnees des noeuds du maillage
%
% function h = POLYFE(n,p,'elem')
% n : matrice N-by-2 
% chaque ligne correspondant aux noeuds d'un element

if nargin==0
    h=struct();
    param=struct();
    domain=[];
    hp = RANDPOLY('fe',param,domain);
    h = class(h,'POLYFE',hp);
    superiorto('RANDPOLY');
else

    if nargin<3
        type='node';
    elseif ~strcmp(type,'elem') && ~strcmp(type,'node')
        warning('le troisieme argument ne veut rien dire')
    end

    h=struct();
    if isa(n,'POLYFE')
        n=getparam(n,'I');
    end
    switch type
    case 'node'        
        if length(n)==1
            param.n=n;
            param.I = [[0:1/n:1-1/n]',[1/n:1/n:1]'];
        else
            if n(1)~= 0 || n(end)~=1
                warning('interval must be [0,1] for polyfe')
                end
                param.I=n(:);
                param.I=[param.I(1:end-1),param.I(2:end)];

            end

        case 'elem'
            param.I=n;

        end
        param.n=size(param.I,1);       
        param.dx = param.I(:,2)-param.I(:,1);


        if nargin==1 || isempty(p)
            p=1;
        end
        param.p=p;
        domain = [0,1];
        hp = RANDPOLY('fe',param,domain);
        h = class(h,'POLYFE',hp);
        superiorto('RANDPOLY');

    end
