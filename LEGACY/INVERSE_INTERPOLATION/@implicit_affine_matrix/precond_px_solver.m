function [u,flag,relres,iter,resvec]=precond_px_solver(P,A,b,k)
%  [u,flag,relres,iter,resvec]=precond_px_solver(P,A,b,k)
% Compute the solution of P(k)A(k) u = P(k)b using cgs solver


tol = 1e-8;
Atmp = eval_sparse(A,k);
afun = @(v) Atmp*v ;
pfun = @(v) cell2mat(cellfun( @(L,U) U\(L\v) ,P.L,P.U,'UniformOutput',0)') * P.Phi(:,k) ;

[u,flag,relres,iter,resvec] = cgs(afun,b,tol,[],pfun);
