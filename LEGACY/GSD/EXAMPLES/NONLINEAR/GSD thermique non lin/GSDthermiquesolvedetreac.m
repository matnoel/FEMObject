function U=GSDthermiquesolvedetreac(S,tol,ka,kn,kl,varargin)

display_ = ischarin('display',varargin);
if display_
    fprintf('\n RESOLUTION PROB DETERMINISTE\n')
end

thermique_nonlin_forms

m = size(ka,1);
U0 = zeros(getnbddlfree(S),m);
U=U0;
n = getnbddlfree(S);

A0 = a{S}(:,:);
b0 = l{S}(:);
if nbsec==2
    b02 = l2{S}(:);
end

A = cell(m,m);
b = cell(m,1);
for i=1:m
    if nbsec==1
        b{i} =  b0*kl(i) ;
    else
        b{i} =  b0*kl{1}(i) +b02*kl{2}(i);
    end
    for j=1:m
        A{i,j} = A0*ka(i,j);
    end
end
A = MULTIMATRIX(A,[n,n],[m,m]);
A = assembleblock(A);
b = MULTIMATRIX(b,[n,1],[m,1]);
b = assembleblock(b);

inc = 0;
for k=1:100
    res = b - A*U(:);
    restemp = cell(m,1);
    for i=1:m
        restemp{i} = zeros(n,1);
        for j=1:m
            for l=1:m
                restemp{i} = restemp{i} + g{S}(:,U(:,j),U(:,l),U*kn(:,l,j,i));
            end
        end
    end
    restemp = MULTIMATRIX(restemp,[n,1],[m,1]);
    restemp = assembleblock(restemp);
    res = res-restemp;
    
    nres = norm(res)/norm(b);
    
    
    if display_
        fprintf('iteration %d , residual = %.2d\n',k,nres);
    end
    if nres<tol
        break
    end
    
    
    AT = A;
    ATtemp = cell(m,m);
    for i=1:m
        for j=1:m
            ATtemp{i,j} =sparse(n,n);
            for l=1:m
                ATtemp{i,j} = ATtemp{i,j} + g{S}(:,:,U(:,l),U*kn(:,l,j,i))+...
                    2*g{S}(:,U(:,l),:,U*kn(:,l,j,i));
            end
        end
    end
    ATtemp = MULTIMATRIX(ATtemp,[n,n],[m,m]);
    ATtemp = assembleblock(ATtemp);
    AT = A + ATtemp;
    %AT = A + dnapprox{S}(:,:,U,U) ;
    %AT=A;
    dU = reshape(AT\res,n,m);
    U = U+dU;
    
    %norm(U-U0)/norm(U)
    
    U0=U;
end

if nres>=tol
    warning('NON CONVERGENCE');
end
