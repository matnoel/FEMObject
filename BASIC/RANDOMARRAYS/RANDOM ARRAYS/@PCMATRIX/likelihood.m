function L = likelihood(Xc,Xs,N,varargin)

% function L = likelihood(X,Xs)
% X : PCMATRIX de taille n-by-1  (vecteur aleatoires)
% Xs : echantillon : matrice de taille n-by-N (N = nombre d'echantillons)
% evalue la log vraisemblance de Xs vis a vis de la mesure de proba de X
%     L = log(prod_k p_X(Xs(:,k)))
% en fait seulement une approximation de ca
%     L ~= log(prod_k prod_i p_Xi(Xs(i,k))) 
%
% function L = likelihood(X,Xs,m)
% m designe le nombre d'echantillons pour l'estimation des densites de
% proba marginales (Par defaut (10^4))
%
% function L = likelihood(X,Xs,m,varargin)
% varargin : arguments passees à la fonction ksdensity

if nargin<3 || isempty(N)
   N=1e4 ;
end

L = 0;
for i=1:numel(Xc)
Xcs = random(getcompo(Xc,i),1,N); % simulation de Xc pour evaluer sa densite de proba
pXc = ksdensity(Xcs,Xs(i,:),varargin{:});
L = L + sum(log(pXc+eps));
end
%L = L / N ; 
return


