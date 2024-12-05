function [Xc,L] = PCidentification_ML(PC,Xs,varargin)
% function [Xc,L] = PCidentification_ML(PC,Xs)
% identification of a random vector on the polynomial chaos
% PC : POLYCHAOS
% Xs : n-by-N array (N realisations d'un vecteur aleatoire de taille n)
%
% function Xc = PCidentification(PC,Xs,'nbsestim',nbsestim)
% 'nbsestim' : nombre de realisations utilisees pour l'estimation des pdf 
% marginales dans l'estimation de la vraisemblance approchee
%
% ---------------------------------------
% CHOIX DE L'OPTIMISATION
% ----------------------------------------
%
% function Xc = PCidentification(PC,Xs,'rs',nbrandini)
% 'rssamples' : nombre de samples du random search 
%
% ----------------------------------------
% EXPLICATION DE LA METHODE UTILISEE      
% Utilisation du maximum de vraisemblance
% ----------------------------------------
% Le vecteur aleatoire est transforme pour avoir une moyenne nulle et 
% une matrice de correlation identite.
% La recherche des coefficients du chaos (n-by-(P+1) coefficients) 
% se ramene donc a la recherche 
% d'une matrice orthogonale n-by-P notee A : on a donc un probleme d'optimisation 
% sur une variete de Stiefel compacte). 
%
%       max(L(A)) avec A matrice orthogonale
%
%   Voir l'article de Soize, CMAME 2010 
%        doi:10.1016/j.cma.2010.03.013

%% definition de fonctions et de la fonction a minimiser
N = size(Xs,2);
mu = sum(Xs,2)/N;
n = size(Xs,1);
Xtildes = (Xs - repmat(mu,1,N));
[U,S,V] = svd(Xtildes,'econ') ;
s = diag(S);
rep = find(s/max(s)>1e-12);
if length(rep)<n
    fprintf('\n !!! variables inter-dependantes ')
    fprintf('\n --> reduction du nombre de variables aleatoires a identifier')    
    n = length(rep);
    U = U(:,rep);
    S = S(rep,rep);
    V=V(:,rep);

    [Xc,L] = PCidentification_ML(PC,S*V',varargin{:}); 
    Xc = mu + U*Xc;
    return
end

U = U*S/sqrt(N);
Ys = V*sqrt(N);
if isa(PC,'double')
    PC = POLYCHAOS(n,PC);
end

P = length(PC)-1;
%dimvariete = n*P-n*(n+1)/2;

fprintf('\n------------------------------------------------------')
fprintf('\n  IDENTIFICATION VARIABLES ALEATOIRES SUR LE CHAOS ')
fprintf('\n------------------------------------------------------')
fprintf('\nNombres de variables a identifier : %d',n)
fprintf('\nNombres d''echantillons            : %d',N)
fprintf('\ndimension du chaos                : %d',getM(PC))
fprintf('\ndegre du chaos                    : %d',max(getorder(PC)))
fprintf('\nnombre de parametres a identifier : %d',n*P)
fprintf('\nnombre de contraintes             : %d',n*(n+1)/2)
fprintf('\n------------------------------------------------------\n')


%% random search  pour l'initialisation
rssamples = getcharin('rssamples',varargin,1000);
Nbs = getcharin('nbsestim',varargin,1e4);

fprintf('\n-----------------------')
fprintf('\n    RANDOM SEARCH    ')
fprintf('\n-----------------------\n')
fprintf('nombre de samples = %d\n',rssamples)


A = zeros(P,n);
perm = 1:n;
%perm = perms(1:m);
Lmin = inf;
for kkk=1:size(perm,1);
    for i=perm
        fprintf('Dimension %d: ',i)
        B = null(A(:,1:i-1)');
        yfun = @(a) PCMATRIX([0,a(:)'],[1,1],PC);
        myfunvrais = @(a,N) -likelihood(yfun(a),Ys(:,i)',N);    
        Pi = P - i +1 ; 
        lb = -ones(Pi,1);
        ub = ones(Pi,1);
        gene = @() B*normalize(lb+rand(Pi,1).*(ub-lb));  
        A(:,i) = my_random_search(@() gene(),@(x) myfunvrais(x,Nbs),rssamples);  
        fprintf('\n')
    end
    YPC = PCMATRIX([zeros(n,1),A'],[n,1],PC);
    L = -likelihood(YPC,Ys',N);    
    if L<Lmin
        YPCmin = YPC ;
        Lmin= L;
    end
end

Xc = mu + U*YPCmin;

fprintf('\nvraisemblance apres random search  = %f\n',L)

fprintf('\n------------------------------------------------------')
fprintf('\n                FIN DE L''IDENTIFICATION          ')
fprintf('\n------------------------------------------------------\n')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FONCTIONS ADDITIONNELLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,L,L0]=my_random_search(argfun,minfun,n)

L = Inf;

for i=1:n
    pourcentage(i,n)
    X0 = argfun();
    L0(i) = minfun(X0);
    s.fval = L0(i);
    if L0(i) < L
        X = X0;
        L = L0(i);
    end
end

return 

function [X,L,L0]=my_lhsrandom_search(lb,ub,minfun,n,illustrate)

r = lhsdesign(n,length(lb));
lb = repmat(lb(:)',n,1);
ub = repmat(ub(:)',n,1);
r = lb+r.*(ub-lb);

L = Inf;

for i=1:n
    pourcentage(i,n)

    X0 = r(i,:);
    L0(i) = minfun(X0);
    s.fval = L0(i);
    marker = 'kx';
    if L0(i) < L
        X = X0;
        L = L0(i);
        marker ='gx';
    end
    state = 'iter';
    if illustrate
        outfun(X0,s,state,marker,minfun)  ;
    end
end


return

