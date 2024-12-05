function [u,result] = solve_alterne_test(A,b,solver,varargin)
% function [u,result] = solve_alterne_test(A,b,solver,varargin)

dim = getdim(A);
param = getparam(solver);

n = size(b);
repone = find(n==1);
param.updatedim(repone)=[];
bu = b;

if ischarin('initialguess',varargin)
    u = getcharin('initialguess',varargin);
    bu = b - A*u;
else
    u = SEPMATRIX(dim);
end

u0=u;
utilde = SEPMATRIX(dim);
hasconv = 0;
erriter = zeros(1,param.maxiter);
errorder = zeros(1,param.maxorder);

if ~isempty(param.metric)
    Ametric = param.metric;
else
    Ametric = cell(1,dim);
    for j=1:dim
        Ametric{j} = speye(n(j));
    end
    Ametric = SEPMATRIX(Ametric);
end

% if dim==2
%     param.maxorder = min(min(n),param.maxorder);
% end

switch param.errorindicator
    case 'reference'
        try
            normref = normND(double(param.reference),Ametric);
        catch
            normref = normND((param.reference),Ametric);
        end
    case 'residual'
        normref = norm(b);
end

for i=1:param.maxorder
    
    hasconviter = 0;
    
    while ~hasconviter
        
        if param.display
            fprintf('order #%d \n',i)
        end
        
        U0 = seprand(n);
        U0tilde = U0;
        alpha0=0;
        
        U = U0;
        Utilde = U0tilde;
        AU = A*U;
        UAU = Utilde'*A*U;
        Ubu = Utilde'*bu;
        
        for kkk=1:param.maxiter    %%%%  ITERATION  ALTERNEE
            
            for j=1:dim %%%% BOUCLE DIMENSION j
                
                
                UAU.F(:,j)= A.F(:,j);
                M = reduce(UAU,j);
                Ubu.F(:,j) = bu.F(:,j);
                S = reduce(Ubu,j);
                
                U.F{j} = M\S;
                alpha = norm(U.F{j});
                
                U.F{j} = U.F{j}/alpha;
                
                if param.adjoint == 1
                    if ismember(j,param.adjointdim)
                        switch param.adjointtype
                            case 0
                                btilde = Ametric*U;
                                Ubtilde = U'*btilde;
                                Ubtilde.F(:,j) = btilde.F(:,j);
                            case 1
                                btilde = Ametric*(alpha*U + u) ;
                                Ubtilde = U'*btilde;
                                Ubtilde.F(:,j) = btilde.F(:,j);
                            case 2
                                btilde = Ametric*(alpha*U + u) - (A')*utilde ;
                                Ubtilde = U'*btilde;
                                Ubtilde.F(:,j) = btilde.F(:,j);
                                
                            otherwise
                                error('pas prevu')
                        end
                        
                        S = reduce(Ubtilde,j);
                        Utilde.F{j} = (M')\S;
                    else
                        Utilde.F{j} = U.F{j};
                    end
                    alphatilde = norm(Utilde.F{j});
                    Utilde.F{j} = Utilde.F{j}/alphatilde;
                else
                    alphatilde = alpha;
                    Utilde = U;
                end
                
                for ii=1:A.m
                    AU.F{ii,j} = A.F{ii,j}*U.F{j};
                    UAU.F{ii,j} = Utilde.F{j}'*AU.F{ii,j};
                end
                for ii=1:bu.m
                    Ubu.F{ii,j} = Utilde.F{j}'*Ubu.F{ii,j};
                end
                
            end %%%% END BOUCLE DIMENSION j
            
            erriter(kkk)=abs(alpha-alpha0)/(alpha+alpha0);
            result.alpha{i}(kkk)=alpha;
            result.alphatilde{i}(kkk)=alphatilde;
            result.U{i}{kkk}=U;
            result.Utilde{i}{kkk}=Utilde;
            
            if param.display
                fprintf('  iteration #%d - stagnation = %d\n',kkk,erriter(kkk))
            end
            
            if erriter(kkk)<param.itercrit
                hasconviter = 1;
                break
            end
            alpha0=alpha;
            
            U0 = U;
            U0tilde=Utilde;
            
            if param.restartifnotconverged==0
                hasconviter=1;
            end
            
        end  %%%%  END ITERATION ALTERNEE
        
    end
    
    
    u = u + alpha*U;
    if param.adjoint
        utilde = utilde + alphatilde*Utilde;
    end
    
    if param.update>0 && i>1 && mod(i,param.updatestep)==0
        
        for kkk=1:param.update
            if param.display
                fprintf(' update #%d, dim ',kkk)
            end
            u0 = u;
            
            W = gathervectors(u);
            if param.adjoint==0 || ~param.updateadjoint
                Wtilde = W;
            else
                Wtilde = gathervectors(utilde);
            end
            
            
            if dim==2 && param.ortho
                W = myorth(W,1:dim);
                Wtilde = myorth(Wtilde,1:dim);
            end
            
            for jj=1:length(param.updatedim)
                j = param.updatedim(jj);
                if param.display
                    fprintf('#%d ',j)
                end
                % if n(jj)<u.m
                %     break
                % end
                
                W.alpha = ones(1,u.m);
                WAW = Wtilde'*A*W;
                Wb = Wtilde'*b;
                WAW.F(:,j)= A.F(:,j);
                Wb.F(:,j)=b.F(:,j);
                M = assembleblock(WAW,j);
                M = M+speye(size(M,1))*param.updateeps;
                S = assembleblock(Wb,j);
                Wj = reshape(M\S,n(j),u.m);
                W.F{1,j}=Wj;
                
                if param.adjoint && param.updateadjoint
                    u = splitvectors(W);
                    u = normalizefuns(u);
                    Wtilde.alpha=ones(1,u.m);
                    AT=A';
                    btemp = Ametric*u;
                    Wb = W'*btemp;
                    WAW = W'*AT*Wtilde;
                    WAW.F(:,j)= AT.F(:,j);
                    Wb.F(:,j)=btemp.F(:,j);
                    M = assembleblock(WAW,j);
                    M = M+speye(size(M,1))*param.updateeps;
                    S = assembleblock(Wb,j);
                    Wj = reshape(M\S,n(j),u.m);
                    Wtilde.F{1,j}=Wj;
                else
                    Wtilde = W;
                end
                
            end
            
            u = splitvectors(W);
            u = normalizefuns(u);
            utilde = splitvectors(Wtilde);
            utilde = normalizefuns(utilde);
            
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
        bu = b - A*u;
        
    else
        bu = bu - alpha*AU;
    end
    
    if param.alphaupdate
        if param.display
            fprintf(' alphaupdate\n')
        end
        W = gathervectors(u);
        W.alpha = ones(1,u.m);
        if param.adjoint
            Wtilde = gathervectors(utilde);
            Wtilde.alpha = ones(1,u.m);
            WAW = timesblock(Wtilde'*A*W);
            Wb =  timesblock(Wtilde'*b);
        else
            WAW = timesblock(W'*A*W);
            Wb =  timesblock(W'*b);
        end
        W.alpha = (WAW\Wb)';
        u=splitvectors(W);
        bu = b - A*u;
        
    end
    
    switch param.errorindicator
        case 'reference'
            try
                errorder(i)=normND(expand(u)-double(param.reference),Ametric)/normref;
            catch
                errorder(i)=normND((u)-(param.reference),Ametric)/normref;
            end
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


function M = reduce(A,j)

Z = A.F(:,j);
A.F(:,j) = [];
a = A.alpha.*prod(cell2mat(A.F),2)';
M = a(1)*Z{1};
for i=2:A.m
    M = M + a(i)*Z{i};
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
