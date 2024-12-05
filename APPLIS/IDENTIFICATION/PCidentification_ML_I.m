function [Xc] = PCidentification_ML_I(PC,Xs,varargin)
% function [Xc,L] = PCidentification_ML_I(PC,Xs)
% identification of a random vector on the polynomial chaos
% PC : POLYCHAOS
% Xs : n-by-N array (N realisations d'un vecteur aleatoire de taille n)
%
% function Xc = PCidentification_ML_I(PC,Xs,'nbsestim',nbsestim)
% 'nbsestim' : nombre de realisations utilisees pour l'estimation des pdf 
% marginales dans l'estimation de la vraisemblance approchee
%
% ----------------------------------------
% EXPLICATION DE LA METHODE UTILISEE      
% Utilisation du Maximum Vraisemblance 
% avec Independece Hypothesis
% ----------------------------------------
% Le vecteur aleatoire est transforme pour avoir une moyenne nulle et 
% une matrice de correlation correlation identite avec Karhunen-Loeve Emprique.
% La recherche des coefficients du chaos (n-by-(P+1) coefficients) 
% se ramene donc a la recherche par dimension, car ici on suppose 
% que les components d'un vector aleatoire sont independents.
% Pour identification d'un dimension, on utilise PCidentification_ML
%
%
%   Voir l'article de Stefanou et Nouy, IJNME 2009
%        doi:10.1002/nme.2546


%% Empirical Karhunen-Loeve expansion 
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
    V = V(:,rep);
end

U = U*S/sqrt(N);
Ys = V*sqrt(N);

if isa(PC,'double')
    PC = POLYCHAOS(n,PC);
end

P = length(PC)-1;

fprintf('\n------------------------------------------------------')
fprintf('\n  IDENTIFICATION VARIABLES ALEATOIRES SUR LE CHAOS  ')
fprintf('\n  Maximum Vraisemblance INDEPENDECE HYPOTHESIS ')
fprintf('\n------------------------------------------------------')
fprintf('\nNombres de variables a identifier : %d',n)
fprintf('\nNombres d''echantillons            : %d',N)
fprintf('\ndimension du chaos                : %d',getM(PC))
fprintf('\ndegre du chaos                    : %d',max(getorder(PC)))
fprintf('\nnombre de parametres a identifier : %d',n*max(getorder(PC)))
fprintf('\n------------------------------------------------------\n')

%% identification of randopm variables by ML avec Independece
% ici on identifie par dimension, et apres on projecte lambda{i}
% rendu par **PCidentification_ML(max(getorder(PC)),Ys(:,i)',varargin{:});** qui 
% est une expansion de chaos de dim 1 et order p
% dans un chaos de dim n et order p. Cependant, que les
% coefficients correspondant chaque dim i sont utilise, les restes sont tous
% zeros. Donc, on peut voir PC qui est dim n et order p comme un chaos
% de n fois de dim 1 et order p. 
lambda = cell(1,n);
Xc = mu;
for i=1:n
    lambda{i} = PCidentification_ML(max(getorder(PC)),Ys(:,i)',varargin{:});
    lambda{i} = setnumber(lambda{i},i);
    lambda{i} = project(lambda{i},PC);
    Xc = Xc + U(:,i)*lambda{i};
end

fprintf('\n------------------------------------------------------')
fprintf('\n                FIN DE L''IDENTIFICATION          ')
fprintf('\n------------------------------------------------------\n')


end
