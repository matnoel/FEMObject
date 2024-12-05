function [u,flag,relres,iter,resvec] = mldivide(X,b,k)
% [u,flag,relres,iter,resvec] = mldivide(X,b,k)
% Compute X(k)\b using cgs solver

error('ZOB')

tol = 1e-8;
afun = @(v) mtimes_truncation(X,v) ;
mfun = @(v) cell2mat(cellfun( @(L,U) L*(U*v) ,X.L,X.U,'UniformOutput',0)') * X.a ;


[u,flag,relres,iter,resvec] = cgs(afun,b,tol,[],mfun);

end

