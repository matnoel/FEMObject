function [u,result,W,Lambda] = solve_gsd_alterne(A,b,solver,solverlambda,solverlambdaupdate)
% function [u,result,W,Lambda] = solve_gsd_alterne(A,b,solver,solverlambda,solverlambdaupdate)

dim = getdim(A);
param = getparam(solver);

n = size(b);
bu = b;

u = SEPMATRIX(dim);
W = zeros(n(1),0);
Lambda = SEPMATRIX(dim);

erriter = zeros(1,param.maxiter);
errorder = zeros(1,param.maxorder);

switch param.errorindicator
    case 'reference'
        normref = normND(param.reference);
    case 'residual'
        normref = norm(b);
end


for i=1:param.maxorder
    fprintf('order #%d \n',i)
    
    lambda = seprand([1,n(2:end)]);
    lambda.F{1,1}=1;
    alpha0=0;
    
    for kkk=1:param.maxiter    %%%%  ITERATION  ALTERNEE
        
        lAl = lambda'*A*lambda;
        lbu = lambda'*bu;
        M = reduce(lAl,1);
        S = reduce(lbu,1);
        U = M\S;
        alpha = norm(U);
        U = U/alpha;
        
        AU = mtimes(A,U,1);
        UAU = mtimes(U',AU,1);
        Ubu = mtimes(U',bu,1);
        
        if dim>2
            [lambda,flag] = solve(UAU,Ubu,solverlambda);
        else
            M = reduce(UAU,2);
            S = reduce(Ubu,2);
            lambda.F{2} = M\S;
        end
        
        
        erriter(kkk)=abs(alpha-alpha0)/(alpha+alpha0);
        fprintf('  iteration #%d - stagnation = %d\n',kkk,erriter(kkk))
        
        if erriter(kkk)<param.itercrit
            break
        end
        alpha0=alpha;
        U0 = U;
        
    end  %%%%  END ITERATION ALTERNEE
    
    % alpha = norm(lambda);
    % stoalpha= [stoalpha, alpha];
    W = [W,full(U)];
    Ulambda = mtimes(U,lambda,1);
    Ulambda = normalizefuns(Ulambda);
    u = u+Ulambda;
    bu = bu - A*Ulambda;
    
    Lambda0 = Lambda;
    for ii=1:Lambda0.m
        Lambda0.F{ii,1} = [Lambda0.F{ii,1};0];
    end
    uni = zeros(i,1);uni(end)=1;
    lambda0 = mtimes(uni,lambda,1);
    Lambda0 = Lambda0+lambda0;
    Lambda = Lambda0;
    
    if param.update && mod(i,param.updatestep)==0 && i>1
        
        fprintf('order #%d - update\n',i)
        % W = orth(W);
        AW = mtimes(A,W,1);
        WAW = mtimes(W',AW,1);
        Wb = mtimes(W',b,1);
        Lambda0=SEPMATRIX(dim);
        
        if dim>2
            Lambda = solve(WAW,Wb,solverlambdaupdate,'initialguess',Lambda0);
        else
            M = assembleblock(WAW,2);
            S = assembleblock(Wb,2);
            Lambda = SEPMATRIX({eye(size(W,2)),reshape(M\S,n(2),size(W,2))});
            Lambda = splitvectors(Lambda);
            Lambda.alpha=ones(1,size(W,2));
        end
        
        u = mtimes(W,Lambda,1);
        bu = b - A*u;
    end
    
    switch param.errorindicator
        case 'reference'
            errorder(i)=normND(expand(u)-param.reference)/normref;
        case 'residual'
            errorder(i)=norm(bu)/normref;
        otherwise
            errorder(i)= norm(Ulambda)/norm(u);%
            % sqrt(min(stoalpha)^2/sum(stoalpha.^2));
    end
    
    fprintf('order #%d - error = %d \n',i,errorder(i))
    
    if errorder(i)<param.tol
        hasconv = 1;
        break
    end
    
end


result.error = errorder;


function M = reduce(A,j)

Z = A.F(:,j);
A.F(:,j) = [];
a = A.alpha.*prod(cell2mat(A.F),2)';
M = a(1)*Z{1};
for i=2:A.m
    M = M + a(i)*Z{i};
end

function M = assembleblock(A,j)

Z = A.F(:,j);
A.F(:,j) = [];
n = size(Z{1},1);
nb = size(A.F{1},1);
B = cell(A.m,1);
for i=1:A.m
    B{i} = A.F{i,1};
    for k=2:size(A.F,2)
        B{i}=B{i}.*A.F{i,k};
    end
end

if size(Z{1},2)>1
    M = sparse([],[],[],n,n,nnz(Z{1})*nb^2);
    for ii=1:nb
        I = (ii-1)*n+(1:n);
        for jj=1:nb
            J = (jj-1)*n+(1:n);
            MIJ = (A.alpha(1)*B{1}(ii,jj))*Z{1};
            for i=2:A.m
                MIJ = MIJ + (A.alpha(i)*B{i}(ii,jj))*Z{i};
            end
            M(I,J)=MIJ;
        end
    end
else
    M = zeros(nb*n,1);
    for ii=1:nb
        I = (ii-1)*n+(1:n);
        MI = (A.alpha(1)*B{1}(ii))*Z{1};
        for i=2:A.m
            MI = MI + (A.alpha(i)*B{i}(ii))*Z{i};
        end
        M(I)=MI;
    end
end


return


function M = timesblock(A)

i=1;
Mi = A.F{i,1};
for k=2:A.dim
    Mi = Mi.*A.F{i,k};
end
M = Mi*A.alpha(i);

for i=2:A.m
    Mi = A.F{i,1};
    for k=2:A.dim
        Mi = Mi.*A.F{i,k};
    end
    M = M + Mi*A.alpha(i);
end

return

