function LU = initialize(LU,A)
% function LU = initialize(N,A)
% LU : LUSOLVER
% A : matrice
% La matrice est factorisee sous la forme A = L*U
% avec L matrice triangulaire inferieure et U matrice triangulaire superieure
% L et U sont stockees et serviront pour la fonction solve

[LU.L,LU.U] = lu(A);
