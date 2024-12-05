function [u,result] = multisvd(b,varargin)
% function [u,result] =  multisvd(b,varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~(length(varargin)==2 && isa(varargin{1},'struct') && isa(varargin{2},'HSEPSOLVER')) %Premier passage
    PP=PARAMETERS(varargin{:});
    if isparamin(PP,'tree')                % L'utilisateur specifie un arbre
        solver = HSEPSOLVER(varargin{:});
        param  = solver(1);
        b=HSEPMATRIX(b,param.tree);
    else                                   % L'utilisateur n'a pas specifie d'arbre
        solver = HSEPSOLVER('tree',b.tree,varargin{:});
        param  = solver(1);
    end
    
else % Recurence...
    param  = varargin{1};
    solver = varargin{2};
end
Tp=param.tree;
SubNodes=find(getconnect(Tp)==param.node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tb  = b.tree;
dim = b.dim;
n   = size(b);
bu  = b;
u   = HSEPMATRIX(Tb);

erriter  = zeros(1,param.maxiter);
errorder = zeros(1,param.maxorder);

switch param.errorindicator
    case {'residual','reference'}
        normref = norm(b);
end

for i=1:param.maxorder
    if param.display
        fprintf('order #%d \n',i)
    end
    
    
    %%% Initialisation du nouveau rang
    U0 = hseprand(Tb,n);
    alpha0=0;
    U = U0;
    
    % Les scalaires : (a optimiser!)
    UbuF=zeros(bu.m,bu.dim);
    UbuA=zeros(1,bu.m);
    for ii=1:bu.m
        for d=1:bu.dim
            UbuF(ii,d)=fastprodscal(U.F{1,d},bu.F{ii,d});
        end
        UbuA(ii)=U.alpha(1)*bu.alpha(ii);
    end
    
    for kkk=1:param.maxiter    %%%%  ITERATION  ALTERNEE
        
        for j=1:dim %%%% BOUCLE DIMENSION j
            
            %%% REDUCE
            S = reduce(UbuF,UbuA,j,bu.F(:,j));
            
            
            %%% STARTWITH ?
            %             if param.start==1
            %                 subparam.startwith=U.F{1,j};
            %             end
            
            %%% RECURCIV'
            U.F{1,j} = multisvd(S, solver(SubNodes(j)) ,...
                setparam(solver,'node',SubNodes(j))  );
            
            %%% ACTUALISER
            alpha = norm(U.F{1,j});
            U.F{1,j} = U.F{1,j}*(1/alpha);
            for ii=1:bu.m
                UbuF(ii,j) = fastprodscal(U.F{1,j},bu.F{ii,j});
                %                 UbuA(ii)   = U.alpha(1)*bu.alpha(ii);
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
    bu = bu - alpha*U;
    
    if param.update>0 && i>1 && mod(i,param.updatestep)==0
        for kkk=1:param.update
            if param.display
                fprintf('update #%d, dim ',kkk)
            end
            
            
            
            
            
            %%% test sortie...
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
        bu = b - u;
    end
    
    
    if param.alphaupdate && mod(i,param.updatestep)==0
        if param.display;fprintf(' alphaupdate\n');end
        u  = solve_alphaupdate(b,u);
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



function M = reduce(F,A,j,Z)
% Compression de l'info :
%   produit des autres dimension
%   somme sur les rangs
F(:,j)=1;
a = (A).*prod(F,2)';
M = a(1)*Z{1};
for i=2:length(a)
    M = M + a(i)*Z{i};
end
return


