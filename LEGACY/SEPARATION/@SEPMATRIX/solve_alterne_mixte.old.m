function [u,result] = solve_alterne_mixte(A,b,solver,varargin)
% function [u,result] = solve_alterne_mixte(A,b,solver,varargin)

A11 = A{1,1};
A12 = A{1,2};
A21 = A{2,1};
A22 = A{2,2};
b1=b{1};
b2=b{2};
dim = getdim(A11);
param = getparam(solver);

n1 = size(b1);
n2 = size(b2);

repone1 = find(n1==1);
repone2 = find(n2==1);
param.updatedim(repone1)=[];
bu1 = b1;
bu2 = b2;


u1 = SEPMATRIX(dim);
u2 = SEPMATRIX(dim);

u01=u1;
u02=u2;
hasconv = 0;
erriter = zeros(1,param.maxiter);
errorder = zeros(1,param.maxorder);

switch param.errorindicator
    case 'reference'
        normref = sqrt(normND(double(param.reference{1}))^2+...
            normND(double(param.reference{2}))^2);
    case 'residual'
        normref = sqrt(norm(b1)^2+norm(b2)^2);
end

for i=1:param.maxorder
    if param.display
        fprintf('order #%d \n',i)
    end
    
    U01 = seprand(n1);
    U02 = seprand(n2);
    alpha01=0;
    alpha02=0;
    
    U1 = U01;
    U2 = U02;
    AU11 = A11*U1;
    AU12 = A12*U2;
    AU21 = A21*U1;
    AU22 = A22*U2;
    
    UAU11 = U1'*AU11;
    UAU12 = U1'*AU12;
    UAU21 = U2'*AU21;
    UAU22 = U2'*AU22;
    Ubu1 = U1'*bu1;
    Ubu2 = U2'*bu2;
    
    for kkk=1:param.maxiter    %%%%  ITERATION  ALTERNEE
        
        for j=1:dim %%%% BOUCLE DIMENSION j
            
            UAU11.F(:,j)= A11.F(:,j);
            UAU12.F(:,j)= A12.F(:,j);
            UAU21.F(:,j)= A21.F(:,j);
            UAU22.F(:,j)= A22.F(:,j);
            M11 = reduce(UAU11,j);
            M12 = reduce(UAU12,j);
            M21 = reduce(UAU21,j);
            M22 = reduce(UAU22,j);
            Ubu1.F(:,j) = bu1.F(:,j);
            Ubu2.F(:,j) = bu2.F(:,j);
            S1 = reduce(Ubu1,j);
            S2 = reduce(Ubu2,j);
            temp = solvesingular([M11,M12;M21,M22],[S1;S2]);
            U1.F{j} = temp(1:n1(j));
            U2.F{j} = temp(n1(j)+1:end);
            alpha1 = norm(U1.F{j});
            alpha2 = norm(U2.F{j});
            U1.F{j} = U1.F{j}/alpha1;
            U2.F{j} = U2.F{j}/alpha2;
            for ii=1:A11.m
                AU11.F{ii,j} = A11.F{ii,j}*U1.F{j};
                UAU11.F{ii,j} = U1.F{j}'*AU11.F{ii,j};
            end
            for ii=1:A12.m
                AU12.F{ii,j} = A12.F{ii,j}*U2.F{j};
                UAU12.F{ii,j} = U1.F{j}'*AU12.F{ii,j};
            end
            for ii=1:A21.m
                AU21.F{ii,j} = A21.F{ii,j}*U1.F{j};
                UAU21.F{ii,j} = U2.F{j}'*AU21.F{ii,j};
            end
            for ii=1:A22.m
                AU22.F{ii,j} = A22.F{ii,j}*U2.F{j};
                UAU22.F{ii,j} = U2.F{j}'*AU22.F{ii,j};
            end
            for ii=1:bu1.m
                Ubu1.F{ii,j} = U1.F{j}'*Ubu1.F{ii,j};
            end
            for ii=1:bu2.m
                Ubu2.F{ii,j} = U2.F{j}'*Ubu2.F{ii,j};
            end
        end %%%% END BOUCLE DIMENSION j
        
        erriter(kkk)=abs(alpha1-alpha01)/(alpha1+alpha01)+abs(alpha2-alpha02)/(alpha2+alpha02);
        if param.display
            fprintf('  iteration #%d - stagnation = %d\n',kkk,erriter(kkk))
        end
        
        if erriter(kkk)<param.itercrit
            break
        end
        alpha01=alpha1;
        alpha02=alpha2;
        U01 = U1;
        U02 = U2;
        
    end  %%%%  END ITERATION ALTERNEE
    
    u1 = u1 + alpha1*U1;
    u2 = u2 + alpha2*U2;
    
    if param.update>0 && i>1 && mod(i,param.updatestep)==0
        
        for kkk=1:param.update
            if param.display
                fprintf(' update #%d, dim ',kkk)
            end
            u01 = u1;
            u02 = u2;
            W1 = gathervectors(u1);
            W2 = gathervectors(u2);
            
            if param.ortho
                W1 = myorth(W1,param.updatedim);
                W2 = myorth(W2,param.updatedim);
            end
            
            for jj=1:length(param.updatedim)
                j = param.updatedim(jj);
                if param.display
                    fprintf('#%d ',j)
                end
                
                W1.alpha = ones(1,u1.m);
                W2.alpha = ones(1,u2.m);
                
                WAW11 = W1'*A11*W1;
                WAW12 = W1'*A12*W2;
                WAW21 = W2'*A21*W1;
                WAW22 = W2'*A22*W2;
                Wb1 = W1'*b1;
                Wb2 = W2'*b2;
                
                WAW11.F(:,j)= A11.F(:,j);
                WAW12.F(:,j)= A12.F(:,j);
                WAW21.F(:,j)= A21.F(:,j);
                WAW22.F(:,j)= A22.F(:,j);
                Wb1.F(:,j)=b1.F(:,j);
                Wb2.F(:,j)=b2.F(:,j);
                M11 = assembleblock(WAW11,j);
                M12 = assembleblock(WAW12,j);
                M21 = assembleblock(WAW21,j);
                M22 = assembleblock(WAW22,j);
                M = [M11,M12;M21,M22] +speye(size(M11,1)+size(M22,1))*param.updateeps;
                S1 = assembleblock(Wb1,j);
                S2 = assembleblock(Wb2,j);
                S = [S1;S2];
                temp = full(M\S);
                W1j = reshape(temp(1:n1(j)*u1.m),n1(j),u1.m);
                W2j = reshape(temp(n1(j)*u1.m+1:end),n2(j),u2.m);
                W1.F{1,j}=W1j;
                W2.F{1,j}=W2j;
            end
            
            u1 = splitvectors(W1);
            u1 = normalizefuns(u1);
            u2 = splitvectors(W2);
            u2 = normalizefuns(u2);
            
            if kkk<param.update
                errup0= sqrt(min(u1.alpha)^2/sum(u1.alpha.^2)) + ...
                    sqrt(min(u2.alpha)^2/sum(u2.alpha.^2));
                % stagn = abs(errup-errup0)/errup;
                stagn = full(sqrt(norm(u1-u01)^2+norm(u2-u02)^2)/sqrt(norm(u01)^2+norm(u02)^2));
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
        bu1 = b1 - A11*u1 - A12*u2;
        bu2 = b2 - A21*u1 - A22*u2;
        
    else
        bu1 = bu1 - alpha1*AU11 - alpha2*AU12;
        bu2 = bu2 - alpha1*AU21 - alpha2*AU22;
    end
    
    
    switch param.errorindicator
        case 'reference'
            errorder(i)=sqrt(normND(expand(u1)-double(param.reference{1}))^2 + ...
                normND(expand(u2)-double(param.reference{2}))^2)/normref;
        case 'residual'
            errorder(i)=sqrt(norm(bu1)^2+norm(bu2)^2)/normref;
        case 'stagnation'
            errorder(i)=sqrt(norm(u1-u01)^2+norm(u2-u02)^2)/...
                sqrt(norm(u1+u01)^2 +norm(u2+u02)^2);
        otherwise
            errorder(i)= sqrt(min(u1.alpha)^2/sum(u1.alpha.^2))+...
                sqrt(min(u2.alpha)^2/sum(u2.alpha.^2));
    end
    if param.display
        fprintf('    order #%d - error = %d \n',i,errorder(i))
    end
    u01=u1;
    u02=u2;
    
    
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

u = {u1;u2};

function M = reduce(A,j)

Z = A.F(:,j);
A.F(:,j) = [];
a = A.alpha.*prod(cell2mat(A.F),2)';
M = a(1)*Z{1};
for i=2:A.m
    M = M + a(i)*Z{i};
end

return

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
    
    M = spalloc(n*nb,n*nb,min((n*nb)^2,sum(nz)*nb^2));
    
    for i=1:A.m
        I = zeros(max(nz),1,nb,nb);
        J = I;
        V = I;
        for ii=1:nb
            for jj=1:nb
                [It,Jt,Vt] = find((A.alpha(i)*B{i}(ii,jj))*Z{i});
                I(1:length(It),1,ii,jj)=It + (ii-1)*n;
                J(1:length(It),1,ii,jj)=Jt + (jj-1)*n;
                V(1:length(It),1,ii,jj)=Vt;
            end
        end
        rep = find(V(:));
        I = I(rep);
        J = J(rep);
        V = V(rep);
        M = M + sparse(I(:),J(:),V(:),n*nb,n*nb);
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


function W = myorth(W,j)

if size(W.F,1)>1
    error('regrouper les fonctions')
end
if nargin==1
    j = 1:size(W.F,2);
end
for k=1:length(j)
    W.F{j(k)} = mygram(W.F{j(k)});
end

return

function W = mygram(W)

for i=1:size(W,2)
    w = W(:,i);
    w = w/sqrt(w'*w);
    w = w - W(:,1:i-1)*(W(:,1:i-1)'*w);
    n = sqrt(w'*w);
    if n<1e-15
        warning('orthogonalisation critique')
    end
    w = w/n;
    W(:,i) = w;
end

return




%
% function M = assembleblock(A,j)
%
% Z = A.F(:,j);
% A.F(:,j) = [];
% n = size(Z{1},1);
% nb = size(A.F{1},1);
% B = cell(A.m,1);
% for i=1:A.m
% B{i} = A.F{i,1};
% for k=2:size(A.F,2)
% B{i}=B{i}.*A.F{i,k};
% end
% end
%
%
%
% if size(Z{1},2)>1
%
% nz = zeros(1,length(Z));
% for i=1:length(Z)
%    nz(i)=nnz(Z{i});
% end
% %M = spalloc(n*nb,n*nb,sum(nz)*nb^2);
% %M = zeros(n*nb);
%
% I = zeros(max(nz),A.m,nb,nb);
% J = I;
% V = I;
% for ii=1:nb
% for jj=1:nb
% for i=1:A.m
% [It,Jt,Vt] = find((A.alpha(i)*B{i}(ii,jj))*Z{i});
% I(1:length(It),i,ii,jj)=It + (ii-1)*n;
% J(1:length(It),i,ii,jj)=Jt + (jj-1)*n;
% V(1:length(It),i,ii,jj)=Vt;
% end
%
% %I(1:length(It),ii,jj) = It;
% %J(1:length(It),ii,jj) = Jt;
% %V(1:length(It),ii,jj) = Vt;
% end
% end
% rep = find(V(:));
% I = I(rep);
% J = J(rep);
% V = V(rep);
% M = sparse(I(:),J(:),V(:),n*nb,n*nb);
%
% else
% M = zeros(nb*n,1);
% for ii=1:nb
% I = (ii-1)*n+(1:n);
% MI = (A.alpha(1)*B{1}(ii))*Z{1};
% for i=2:A.m
% MI = MI + (A.alpha(i)*B{i}(ii))*Z{i};
% end
% M(I)=MI;
% end
% end
%
%
%
% return
%
