function L = likelihood(Xc,Xs)

% function L = likelihood(X,Xs)
% X : RANDVAR
% Xs : echantillon : matrice de taille N-by-1 (N = nombre d'echantillons)
% evalue la vraisemblance de Xs vis a vis de la mesure de proba de X
%     L = sum_k log(p_X(Xs(k)))


pXc = pdf(Xc,Xs);

L = sum(log(pXc+eps));

