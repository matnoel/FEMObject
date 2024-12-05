function h = POLYFELAGRANGE(n,m,varargin)
% function h = POLYFELAGRANGE(n,m)
% base d'interpolation de lagrange-lobatto par morceaux 
% m : nombre de points d'interpolation par morceau
% n : si scalaire : decoupage uniforme de [0,1] en n elements
%     si vecteur  : coordonnes des noeuds du maillage

% function h = POLYFELAGRANGE(n,m,'node')
% n correspond aux coordonnees des noeuds du maillage
%
% function h = POLYFELAGRANGE(n,m,'elem')
% n : matrice N-by-2 
% chaque ligne correspondant aux noeuds d'un element
%
%
% function h = POLYFELAGRANGE(n,m,'elem',pol)
% pol : RANDPOLY servant a la definition de la grille de Lobatto (POLYLEGENGRE, ... )

if nargin==0
    h=struct();
    param=[];
    domain=[];
    hp = RANDPOLY('felagrange',param,domain);
    h = class(h,'POLYFELAGRANGE',hp);
    superiorto('RANDPOLY');
else


    h=struct();

    type = getclassin('char',varargin,'node');
    pol = POLYLEGENDRE();

    fe = POLYFE(n,0,type);
    param = getparam(fe);

    g = calc_gausslobattopoints(pol,m);

    param.points=zeros(m*param.n-(param.n-1),1);
    param.weights=zeros(1,m*param.n-(param.n-1));
    for i=1:param.n     
        x1 = param.I(i,1);
        x2 = param.I(i,2);
        dx=x2-x1;
        param.subpoints{i} = transfer(RVUNIFORM(-1,1),RVUNIFORM(x1,x2),g.coord);
        param.reppoints{i} = (i-1)*(m-1)+(1:m);
        param.points(param.reppoints{i}) = param.subpoints{i}; 
        param.subweights{i} = dx*g.w; 
        param.weights(param.reppoints{i}) = ...
            param.weights(param.reppoints{i})+ param.subweights{i}; 
    end

    param.m=m;
    param.p=m-1;
    domain = [param.points(1),param.points(end)];
    hp = RANDPOLY('felagrange',param,domain);
    h = class(h,'POLYFELAGRANGE',hp);
    superiorto('RANDPOLY');

end
