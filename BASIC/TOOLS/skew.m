function S = skew(n,a)
% function S = skew(n,a)
% cree une matrice antisymetrique de taille n-by-n
% a partir d'un vecteur a de taille n(n-1)/2

if length(a)~=(n*(n-1)/2)
    error('le vecteur doit etre de taille %d',(n*(n-1)/2))
end
rep = reshape(1:n^2,n,n);
repL = nonzeros(tril(rep,-1));

S = zeros(n,n);
S(repL(:))=a;
S = S-S';

