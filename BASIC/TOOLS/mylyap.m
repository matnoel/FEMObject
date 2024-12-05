function X = mylyap(A,B,Q)
% function X = mylyap(A,B,Q)
% Solve AX + XB + Q = 0

n = size(A,1);
I = speye(n);

if n < 100
    s = kron(I,A);
    v = kron(B,I);
    s = s + v;
    
    Qt = Q';
    QQ = Qt(:);
    
    %[L,U] = ilu(s);
    %X = qmr(s,-QQ,5e-2,min(500,size(s,1)));
    X = s\-QQ;
    X = reshape(X,n,n);
    X = X';
else
    D = SEPMATRIX({I,A;B',I});
    [U,S,V] = svd(full(-Q));
    S = diag(S)';
    E = SEPMATRIX({U,V},S);
    E = splitvectors(E);
    
    PGD = SEPSOLVER(2,'maxorder',500,...
        'display',1,'errorindicator','residual','tol',1e-3);
    [X,r] = solve_alterne_V2(D'*D,D'*E,PGD);
    fprintf('PGD converged in %d iterations with a relative residual of %d\n',X.m,r.error(X.m))
    X = expand(X);
end
