function [u,result] = multisvd(b,varargin)
% function [u,result] =  multisvd(b,varargin)

if any(cellfun(@issparse,b.F(:))) % avoid bad surprises...
    b=full(b);
end

dim = getdim(b);
if nargin==3 && isa(varargin{1},'struct') && isa(varargin{2},'HSEPSOLVER') % Cas hierarchique
    param=varargin{1};
else
    solver = SEPSOLVER(dim,varargin{:});
    param = getparam(solver);
end


n = size(b);
bu = b;
u = SEPMATRIX(dim);

erriter = zeros(1,param.maxiter);
errorder = zeros(1,param.maxorder);

switch param.errorindicator
    case {'residual','reference'}
        normref = norm(b);
end


for i=1:param.maxorder
    if param.display
        fprintf('order #%d \n',i)
    end
    
    U0 = seprand(n);
    U0 = normalizefuns(U0);
    alpha0=1;
    U = U0;
    Ubu = U'*bu;
    
    for kkk=1:param.maxiter    %%%%  ITERATION  ALTERNEE
        for j=1:dim %%%% BOUCLE DIMENSION j
            
            Ubu.F(:,j) = bu.F(:,j);
            S = reduce(Ubu,j);
            U.F{j} = S;
            alpha = norm(U.F{j});
            U.F{j} = U.F{j}/alpha;
            global mexcompiled;
            if ~mexcompiled
                for ii=1:bu.m
                    Ubu.F{ii,j} = U.F{j}'*Ubu.F{ii,j};
                end
            else
                Ubu.F(:,j) = multiplyF({U.F{j}'},Ubu.F(:,j));
            end
            
        end %%%% END BOUCLE DIMENSION j
        
        erriter(kkk)=abs(alpha-alpha0)/(alpha+alpha0);
%         erriter(kkk)= norm(alpha*U-alpha0*U0)/alpha;
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
    
    if param.update>0 && i>1 && mod(i,param.updatestep)==0
        for kkk=1:param.update
            if param.display
                fprintf('update #%d, dim ',kkk)
            end
            u0=u;
            W = gathervectors(u);
            if dim==2
                W = orth(W);
            end
            for jj=1:length(param.updatedim)
                j = param.updatedim(jj);
                if param.display
                    fprintf('#%d ',j)
                end
                W.alpha = ones(1,i);
                WW = W'*W;
                Wb = W'*b;
                WW.F(:,j)= {ones(i,i)};
                Wb.F(:,j)= b.F(:,j);
                M = timesblock(WW);
                M = M +1e-14*eye(size(M,1));
                S = reshape(assembleblock(Wb,j),n(j),i);
                Wj = (M\S')';
                W.alpha = ones(1,i);
                W.F{1,j}=Wj;
            end
            u = splitvectors(W);
            u = normalizefuns(u);
            if kkk<param.update
                errup0= sqrt(min(u.alpha)^2/sum(u.alpha.^2));
                %stagn = abs(errup-errup0)/errup;
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
        bu = b - u;
    else
        bu = bu - alpha*U;
    end
    
    
    if param.alphaupdate
        W = gathervectors(u);
        W.alpha = ones(1,i);
        WW = timesblock(W'*W);
        Wb = timesblock(W'*b);
        W.alpha = (WW\Wb)';
        u=splitvectors(W);
        bu = b - u;
    end
    
    switch param.errorindicator
        case {'residual','reference'}
            errorder(i)=norm(bu)/normref;
        otherwise
            errorder(i)= sqrt(min(u.alpha)^2/sum(u.alpha.^2));
    end
    
    if param.display
        fprintf('  order #%d - error = %d \n',i,errorder(i))
    end
    
    if errorder(i)<param.tol
        hasconv = 1;
        break
    end
    
end


result.error = errorder;

%
% function M = reduce(A,j)
%
% Z = A.F(:,j);
% A.F(:,j) = [];
% a = A.alpha.*prod(cell2mat(A.F),2)';
% M = a(1)*Z{1};
% for i=2:A.m
% M = M + a(i)*Z{i};
% end
%








