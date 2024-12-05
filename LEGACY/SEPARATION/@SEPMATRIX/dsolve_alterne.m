function [u,result] = dsolve_alterne(B,A,b,N,solver,varargin)
% function [u,result] = dsolve_alterne(B,A,b,N,solver,varargin)

dim = getdim(A);
param = getparam(solver);

n = size(b);
nx = n(1);
nt = n(2);
Mt = getMmatrix(N);
Dt = getDmatrix(N,'basic');
Dtdiff = getDmatrix(N);

repone = find(n==1);
param.updatedim(repone)=[];
bu = b;

if ischarin('initialguess',varargin)
    error('')
    %   u = getcharin('initialguess',varargin);
    %   bu = b - A*u;
else
    u = SEPMATRIX(dim);
end

u0=u;
utilde = SEPMATRIX(dim);
hasconv = 0;
erriter = zeros(1,param.maxiter);
errorder = zeros(1,param.maxorder);

% if dim==2
%     param.maxorder = min(min(n),param.maxorder);
% end

switch param.errorindicator
    case 'reference'
        normref = normND(param.reference);
    case 'residual'
        normref = norm(b);
end

for i=1:param.maxorder
    if param.display
        fprintf('order #%d \n',i)
    end
    
    U0 = seprand(n);
    alpha0=0;
    
    U = U0;
    Utilde = U;
    BU = B*U;
    AU = A*U;
    UAU = Utilde'*AU;
    UBU = Utilde'*BU;
    Ubu = Utilde'*bu;
    
    for kkk=1:param.maxiter    %%%%  ITERATION  ALTERNEE
        
        UAU.F(:,1) = A.F(:,1);
        UBU.F(:,1) = B.F(:,1);
        Ubu.F(:,1) = bu.F(:,1);
        MB = reduce(UBU,1);
        MA = reduce(UAU,1);
        S = reduce(Ubu,1);
        U.F{1} = full(double(dsolve(N,S,MB,MA)));
        U.F{1} = U.F{1}/norm(U.F{1});
        if param.adjoint == 1
            error('')
        else
            Utilde = U;
        end
        for ii=1:A.m
            AU.F{ii,1} = A.F{ii,1}*U.F{1};
            UAU.F{ii,1} = sum(sum(Utilde.F{1}.*(AU.F{ii,1}*Mt)));
        end
        for ii=1:B.m
            BU.F{ii,1} = B.F{ii,1}*U.F{1};
            UBU.F{ii,1} = sum(sum(Utilde.F{1}.*(BU.F{ii,1}*Dt')));
        end
        for ii=1:bu.m
            Ubu.F{ii,1} = sum(sum(Utilde.F{1}.*(Ubu.F{ii,1}*Mt)));
        end
        
        for j=2:dim %%%% BOUCLE DIMENSION j
            
            UAU.F(:,j)= A.F(:,j);
            M = reduce(UAU,j);
            Ubu.F(:,j) = bu.F(:,j);
            S = reduce(Ubu,j);
            U.F{j} = M\S;
            alpha = norm(U.F{j});
            U.F{j} = U.F{j}/alpha;
            
            if param.adjoint == 1
                error('')
            else
                Utilde = U;
            end
            
            for ii=1:A.m
                AU.F{ii,j} = A.F{ii,j}*U.F{j};
                UAU.F{ii,j} = Utilde.F{j}'*AU.F{ii,j};
            end
            for ii=1:B.m
                BU.F{ii,j} = B.F{ii,j}*U.F{j};
                UBU.F{ii,j} = Utilde.F{j}'*BU.F{ii,j};
            end
            for ii=1:bu.m
                Ubu.F{ii,j} = Utilde.F{j}'*Ubu.F{ii,j};
            end
            
        end %%%% END BOUCLE DIMENSION j
        
        erriter(kkk)=abs(alpha-alpha0)/(alpha+alpha0);
        if param.display
            fprintf('  iteration #%d - stagnation = %d\n',kkk,erriter(kkk))
        end
        
        if erriter(kkk)<param.itercrit
            break
        end
        alpha0=alpha;
        U0 = U;
        
    end  %%%%  END ITERATION ALTERNEE
    
    u = u + alpha*U;
    if param.adjoint
        utilde = utilde + Utilde;
    end
    
    if param.update>0 && i>1 && mod(i,param.updatestep)==0
        
        for kkk=1:param.update
            if param.display
                fprintf(' update #%d, dim ',kkk)
            end
            u0 = u;
            
            Wxt = u.F(:,1);
            W = u;
            W.alpha = ones(1,u.m);
            W.F(:,1)={1};
            W = gathervectors(W);
            W.F{1}=1;
            WAW = W'*A*W;
            WBW = W'*B*W;
            Wb = W'*b;
            WAW.F(:,1)={zeros(i,i)};
            WBW.F(:,1)={zeros(i,i)};
            Wb.F(:,1)={zeros(i,1)};
            
            for ii=1:A.m
                for iii=1:i
                    for jjj=1:i
                        WAW.F{ii,1}(iii,jjj) = sum(sum(Wxt{iii}.*((A.F{ii,1}*Wxt{jjj})*Mt)));
                    end
                end
            end
            for ii=1:B.m
                for iii=1:i
                    for jjj=1:i
                        WBW.F{ii,1}(iii,jjj) = sum(sum(Wxt{iii}.*((B.F{ii,1}*Wxt{jjj})*Dt')));
                    end
                end
            end
            for ii=1:b.m
                for iii=1:i
                    Wb.F{ii,1}(iii,1) = sum(sum(Wxt{iii}.*(b.F{ii,1}*Mt)));
                end
            end
            
            for jj=1:length(param.updatedim)
                W.alpha = ones(1,u.m);
                j = param.updatedim(jj);
                if param.display
                    fprintf('#%d ',j)
                end
                
                WAW.F(:,j)= A.F(:,j);
                WBW.F(:,j) = B.F(:,j);
                Wb.F(:,j)=b.F(:,j);
                M = assembleblock(WAW,j)+assembleblock(WBW,j);
                M = M+speye(size(M,1))*1e-14;
                S = assembleblock(Wb,j);
                Wj = reshape(M\S,n(j),u.m);
                W.alpha = ones(1,u.m);
                W.F{1,j}=Wj;
                for ii=1:A.m
                    WAW.F{ii,j} = Wj'*A.F{ii,j}*Wj;
                end
                for ii=1:B.m
                    WBW.F{ii,j} = Wj'*B.F{ii,j}*Wj;
                end
                for ii=1:b.m
                    Wb.F{ii,j} = Wj'*b.F{ii,j};
                end
                
                % if param.adjoint
                % u = splitvectors(W);
                % u = normalizefuns(u);
                % AT=A';
                % Wb = W'*u;
                % WAW = W'*AT*Wtilde;
                % WAW.F(:,j)= AT.F(:,j);
                % Wb.F(:,j)=u.F(:,j);
                % M = assembleblock(WAW,j);
                % S = assembleblock(Wb,j);
                % Wj = reshape(M\S,n(j),u.m);
                % Wtilde.alpha = ones(1,u.m);
                % Wtilde.F{1,j}=Wj;
                % end
            end
            u.alpha = ones(1,u.m);
            u.F(:,1) = Wxt;
            W.F{1} = eye(u.m);
            W = splitvectors(W);
            u.F(:,2:end) = W.F(:,2:end);
            u = normalizefuns(u);
            
            if kkk<param.update
                errup0= sqrt(min(u.alpha)^2/sum(u.alpha.^2));
                % stagn = abs(errup-errup0)/errup;
                stagn = full(norm(u-u0)/norm(u0));
                if param.display
                    fprintf('- stagnation = %d\n',stagn)
                end
                if stagn<errup0*param.itercritupdate
                    break
                end
            else
                if param.display
                    fprintf('\n')
                end
            end
            
        end
        
        bu = b - A*u - B*mtimes(u,Dtdiff',1);
        
    else
        
        bu = bu - alpha*mtimes(BU,Dtdiff',1) - alpha*AU;
    end
    
    
    switch param.errorindicator
        case 'reference'
            errorder(i)=normND(expand(u)-param.reference)/normref;
        case 'residual'
            errorder(i)=norm(bu)/normref;
        case 'stagnation'
            errorder(i)=norm(u-u0)/norm(u+u0);
        otherwise
            errorder(i)= sqrt(min(u.alpha)^2/sum(u.alpha.^2));
    end
    if param.display
        fprintf('    order #%d - error = %d \n',i,errorder(i))
    end
    u0=u;
    
    if errorder(i)<param.tol
        hasconv = 1;
        break
    end
    
end
if nargout==1
    if errorder(i)>=param.tol
        fprintf('    no convergence - ')
    else
        fprintf('    convergence - ')
    end
    fprintf('order #%d - error = %d\n',i,errorder(i))
end
result.error = errorder;

function M = orthomat(W)



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
    
    nz = zeros(1,length(Z));
    for i=1:length(Z)
        nz(i)=nnz(Z{i});
    end
    % M = spalloc(n*nb,n*nb,sum(nz)*nb^2);
    % M = zeros(n*nb);
    
    I = zeros(max(nz),A.m,nb,nb);
    J = I;
    V = I;
    for ii=1:nb
        for jj=1:nb
            for i=1:A.m
                [It,Jt,Vt] = find((A.alpha(i)*B{i}(ii,jj))*Z{i});
                I(1:length(It),i,ii,jj)=It + (ii-1)*n;
                J(1:length(It),i,ii,jj)=Jt + (jj-1)*n;
                V(1:length(It),i,ii,jj)=Vt;
            end
            
            % I(1:length(It),ii,jj) = It;
            % J(1:length(It),ii,jj) = Jt;
            % V(1:length(It),ii,jj) = Vt;
        end
    end
    rep = find(V(:));
    I = I(rep);
    J = J(rep);
    V = V(rep);
    M = sparse(I(:),J(:),V(:),n*nb,n*nb);
    
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


