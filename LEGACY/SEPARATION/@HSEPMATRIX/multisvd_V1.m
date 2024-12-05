function [u,result] = multisvd(b,varargin)
% function [u,result] =  multisvd(b,varargin)

solver = HSEPSOLVER(varargin{:});
param  = getparam(solver);

% Initialisation de param si RACINE
root=0;
if param.node==0
    param=initparam(param);
    root=1;
end


%%% C'est parti :
dim = getdim(b);
bu = b;
Tb = gettree(b);
u = HSEPMATRIX(Tb);
n = size(b);
% Initialisation de U : U0
U0=cell(1,param.maxorder);
if param.start && root==0
    for i=1:min(param.startwith.m,param.maxorder)
        U0{i}= truncate(param.startwith,i);
        %         U0{i}.alpha=1;
    end
    for i=min(param.startwith.m,param.maxorder)+1:param.maxorder
        U0{i}= hseprand(Tb,n);                                    %%%%
    end
    
else
    for i=1:param.maxorder
        U0{i}= hseprand(Tb,n);                                    %%%%
    end
end



erriter  = zeros(1,param.maxiter);
errorder = zeros(1,param.maxorder);

switch param.errorindicator
    case {'residual','reference'}
        normref = norm(b);
end
%%%
normref = norm(b);
%%%

cvaltern = cell(param.maxorder,param.maxiter,3);

for i=1:param.maxorder
    if 0%root%param.display
        fprintf('order #%d for node #%d \n',i,param.node)
    end
    
    %%% Initialisation du nouveau rang U
    alpha0=0;                                        %%%%
    U = U0{i};                                       %%%%
    
    % Les scalaires :
    UbuF=zeros(bu.m,bu.dim);
    UbuA=zeros(1,bu.m);
    for ii=1:bu.m
        for d=1:bu.dim
            UbuF(ii,d)=fastprodscal(U.F{1,d},bu.F{ii,d});
        end
        UbuA(ii)=U.alpha(1)*bu.alpha(ii);
    end
    
    for kkk=1:param.maxiter    %%%%  ITERATION  ALTERNEE
        %         alpha=U.alpha(1);
        UbuO = diag(UbuA)*UbuF;
        Uold = U;
        for j=1:dim            %%%%  BOUCLE DIMENSION j
            %%% REDUCE
            S = reduce(UbuF,UbuA,j,bu.F(:,j));
            %%% SUBPARAM
            subparam = calc_subparam(param,j,kkk,i);
            if param.start==1
                subparam.startwith=U.F{1,j};% truncate(U.F{1,j},1:1);
            end
            total = [fieldnames(subparam)';struct2cell(subparam)'];
            %%% RECURCIV'
            U.F{1,j} = multisvd(S,total{:});
            
            %%% ACTUALISER
            alpha = norm(U.F{1,j});
            
            %             if isnan(alpha)
            %                 alpha=alpha
            %             end
            
            U.F{1,j} = U.F{1,j}*(1/alpha);
            for ii=1:bu.m
                UbuF(ii,j) = fastprodscal(U.F{1,j},bu.F{ii,j});
                UbuA(ii)   = bu.alpha(ii)*U.alpha(1);
            end
            
        end %%%% END BOUCLE DIMENSION j
        Unew = U;
        UbuN = diag(UbuA)*UbuF;
        
        %%% STOQUER POUR OBSERVATION
        
        %         erriter(kkk)=abs(alpha-alpha0)/(alpha+alpha0);
        %         erriter(kkk)=norm(U1new-U1old)/(alpha);
        erriter(kkk)=norm(UbuO-UbuN)/norm(UbuN);
        
        
        if root || param.display
            cvaltern{i,kkk,1}=erriter(kkk);
            cvaltern{i,kkk,2}=Unew;% (norm(Uold-Unew)/norm(Unew));
            cvaltern{i,kkk,3}=UbuN;% abs(alpha-alpha0)/(alpha);
            fprintf(' iteration #%d - stagnation = %d\n'  ,kkk,cvaltern{i,kkk,1})
            %             fprintf(' iteration #%d - DELTA TETA = %d\n'  ,kkk,cvaltern(i,kkk,2))
            %             fprintf(' iteration #%d - ALPHA STAG = %d\n\n',kkk,cvaltern(i,kkk,3))
            
        end
        
        if erriter(kkk)<param.itercrit
            %             break
        end
        alpha0=alpha;
        
        %         U0 = U;                                  %%%%%%%%%%%%%%%%%%%%
        
    end  %%%%  END ITERATION ALTERNEE
    
    % Remplir le nouveau rang
    u  = u  + alpha*U;
    bu = bu - alpha*U;
    
    switch param.errorindicator
        case {'residual','reference'}
            errorder(i)=norm(bu)/normref;
        otherwise
            errorder(i)= sqrt(min(u.alpha)^2/sum(u.alpha.^2));
    end
    
    if root || param.display
        fprintf('               order #%d - error = %d \n\n',i,errorder(i))
    end
    
    if errorder(i)<param.tol
        hasconv = 1;
        break
    end
end
result.error = errorder;
result.taux_compression  = nbvect(u)/nbvect(b);
result.nbiteraltern = kkk;
result.cvaltern = cvaltern;




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




function param=calc_subparam(param,j,K,I)
%%% trouver le nouveau noeud et ecrire le bon varargin!
% global NIANIA
subnodes = find(get_connect(param.tree)==param.node);
subnode  = subnodes(j); % RESTE A VERIFIER...
param.node     = subnode;
param.maxorder = param.rank(subnode);
param.depth    = param.depth+1;
if param.dyntol
    param.tol = 3^(-(K)); % NIANIA
end
if param.dynitercrit
    param.itercrit = 5e-2;
    param.maxiter  = 10;
end
% param.start = 0;
return



% Fonction appelee par la RACINE !
function param=initparam(param)

T=param.tree;
c=get_connect(T);
param.node = find(c==0);
% Petite vérification : r=1 pour les feuilles SM d'une seule variable :
v=get_var(T);
r=getrank(param.tree);
for i=1:length(c)
    if isRleaf(T,i) || length(v{i})==1
        r(i)=1;
    end
end
param.rank     = r;
param.maxorder = param.rank(param.node);





