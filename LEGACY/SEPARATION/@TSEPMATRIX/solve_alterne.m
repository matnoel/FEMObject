function [u,result] = solve_alterne(A,b,solver,varargin)
% function [u,result] = solve_alterne(A,b,solver,varargin)

global mexcompiled

dim = getdim(A);
param = getparam(solver);

if param.residual
    solver = setparam(solver,'residual',false);
    [u,result] = solve_alterne(A'*A,A'*b,solver,varargin{:});
    return
end

if isa(b,'SEPMATRIX')
    b=TSEPMATRIX(b);
end
n = size(b);
bu = b;

if ischarin('initialguess',varargin)
    u = getcharin('initialguess',varargin);
    bu = orth(b - orth(A*u));
else
    if param.alphaupdate
        u = TSEPMATRIX(dim);
    else
        u=SEPMATRIX(dim);
    end
end

u0=u;

erriter = zeros(1,param.maxiter);
errorder = zeros(1,param.maxorder);

if ~isempty(param.metric)
    if isa(param.metric,'SEPMATRIX')
        Ametric = TSEPMATRIX(param.metric);
    else
        Ametric = param.metric;
    end
else
    Ametric = cell(1,dim);
    for j=1:dim
        Ametric{j} = speye(n(j));
    end
    Ametric = TSEPMATRIX(Ametric);
end

switch param.errorindicator
    case 'reference'
        if isempty(param.reference)
            warning('no reference: error indicator switched to none')
            param.errorindicator = 'none';
        elseif isa(param.reference,'double')
            normref = normND(param.reference);
        else
            if isa(param.reference,'SEPMATRIX')
                normref = norm(TSEPMATRIX(param.reference));
            else
                normref = norm(param.reference);
            end
        end
    case 'residual'
        normref = norm(b);
end

for i=1:param.maxorder
    if param.display
        fprintf('order #%d \n',i)
    end
    switch param.inittype
        case 'rand'
            U0 = tseprand(n);
            U0.alpha(1) = 1;
        case 'one'
            U0 = tsepone(n);
    end
    
    alpha0=norm(U0);
    U = U0;
    AU = A*U;
    UAU=U'*AU;
    Ubu=U'*bu;
    
    for kkk=1:param.maxiter    %%%%  ITERATION  ALTERNEE
        for j=1:dim %%%% BOUCLE DIMENSION j
            UAU.F{j}(:)=A.F{j}(:);
            M=reduce(UAU,j);
            M=reshape(M,n(j),n(j));
            
            Ubu.F{j}(:)=bu.F{j}(:);
            S=reduce(Ubu,j);
            
            U.F{j}{1}=M\S;
            alpha=norm(U.F{j}{1});
            U.F{j}{1}=U.F{j}{1}/alpha;
            if ~mexcompiled
                AU.F{j}=cellfun(@(x) x*U.F{j}{1},A.F{j},'UniformOutput',false);
                UAU.F{j}=cellfun(@(x) U.F{j}{1}'*x,AU.F{j},'UniformOutput',false);
                Ubu.F{j}=cellfun(@(x) U.F{j}{1}'*x,bu.F{j},'UniformOutput',false);
            else
                AU.F{j} = multiplyF(A.F{j},{U.F{j}{1}});
                UAU.F{j} = multiplyF({U.F{j}{1}'},AU.F{j});
                Ubu.F{j} = multiplyF({U.F{j}{1}'},bu.F{j});
            end
        end %%%% END BOUCLE DIMENSION j
        
        erriter(kkk)=abs(alpha-alpha0)/(alpha+alpha0);
        result.alpha{i}(kkk)=alpha;
        
        if param.display
            fprintf('  iteration #%d - stagnation = %d\n',kkk,erriter(kkk))
        end
        
        if erriter(kkk)<param.itercrit
            break
        end
        alpha0=alpha;
        U0 = U;
    end  %%%%  END ITERATION ALTERNEE
    
    
    if ~param.alphaupdate
        u=u+alpha*SEPMATRIX(U);
        bu = orth(bu-alpha*AU);
    else
        u=orth(u+alpha*U); % on peut faire une orthogonalisation iterative
        if (i>1) && mod(i,param.updatestep)==0
            if param.display;fprintf(' alphaupdate\n');end
            u=solve_alphaupdate(A,b,u,[],param);
            bu=orth(b-orth(A*u));
        else
            bu = orth(bu-alpha*AU);
        end
    end
    
    if param.storeiter==1
        result.u{i}=u;
    end
    
    switch param.errorindicator
        case 'reference'
            if isa(param.reference,'double')
                errorder(i)= normND(expand(u)-param.reference)/normref;
            else
                errorder(i)= norm(u-TSEPMATRIX(param.reference))/normref;
            end
        case 'residual'
            errorder(i)=norm(bu)/normref;
        case 'stagnation'
            errorder(i)=norm(u-u0)/norm(u+u0);
        otherwise
            if isa(u,'SEPMATRIX')
                errorder(i)=sqrt(min(abs(u.alpha))^2/sum(u.alpha.^2));
            else
                % temp=double(u.alpha);
                % errorder(i)=abs(temp(end))/norm(temp(:));
                
                % errorder(i)=alpha/norm(u.alpha);
                
                a=tenzeros(size(u.alpha));
                switch dim
                    case 2
                        if u.m(1)~=1
                            a(end,:)=u.alpha(end,:);
                        end
                        if u.m(2)~=1
                            a(:,end)=u.alpha(:,end);
                        end
                        errorder(i)=norm(a)/norm(u.alpha);
                    case 3
                        if u.m(1)~=1
                            a(end,:,:)=u.alpha(end,:,:);
                        end
                        if u.m(2)~=1
                            a(:,end,:)=u.alpha(:,end,:);
                        end
                        if u.m(3)~=1
                            a(:,:,end)=u.alpha(:,:,end);
                        end
                        errorder(i)=norm(a)/norm(u.alpha);
                    otherwise
                        error('not implemented');
                end
            end
    end
    if param.display
        fprintf('    order #%d - error = %d \n',i,errorder(i))
    end
    
    u0=u;
    
    if errorder(i)<param.tol
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
