function [x,flag] = solve(C,A,b)
% function [x,flag] = solve(C,A,b)
% resolution de Ax=b avec factorisation de Cholesky
% C : CHOLESKYSOLVER
% A : matrice
% b : vecteur
% function [x,flag] = solve(C,b)
% resolution de U'*Ux=b apres calcul de la matrice U via l'appel de la
% fonction initialize du solveur C

if nargin==3
    C.U = chol(A);
elseif nargin==2 && isempty(C.U)
    error('rentrer une matrice ou initialiser le solveur')
end

x = C.U\(C.U'\b);
flag = 0;
