function [x,flag] = solve(LU,A,b)
% function [x,flag] = solve(LU,A,b)
% resolution de Ax=b avec factorisation LU
% LU : LUSOLVER
% A : matrice
% b : vecteur
% function [x,flag] = solve(LU,b)
% resolution de L*Ux=b apres calcul des matrices L et U via l'appel de la
% fonction initialize du solveur LU

if nargin==3
    [LU.L,LU.U] = lu(A);
elseif nargin==2 && (isempty(LU.U) || isempty(LU.L))
    error('rentrer une matrice ou initialiser le solveur')
end

x = LU.U\(LU.L\b);
flag = 0;
