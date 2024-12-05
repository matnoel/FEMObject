function C = initialize(C,A)
% function C = initialize(C,A)
% C : CHOLESKYSOLVER
% A : matrice
% La matrice est factorisee sous la forme A = U'*U
% avec U matrice triangulaire superieure
% U est stockee et servira pour la fonction solve

C.U = chol(A);
