function [um,result] = multisvd_constrained(u,a,b,epsilon,varargin)
% function [um,result] =  multisvd_constrained(u,a,b,epsilon,varargin)
% u : ND - array
% a : borne min (double ou ND array), eventuellement []
% b : borne max (double ou ND array), eventuellement []
% epsilon : parametre de penalisation

N = NEWTONSOLVER('type','full','display',false,'tol',1e-5,'maxiter',100);

if isempty(a) && isempty(b)
    fun1 = @(x) zeros(size(x));
    funtang1 = @(x) zeros(size(x));
elseif ~isempty(a) && ~isempty(b)
    fun1 = @(x) epsilon*(max(x-b,0) - max(a-x,0));
    funtang1 = @(x) epsilon*(double(x>b) + double(x<a));
elseif isempty(b)
    fun1 = @(x) epsilon*(-max(a-x,0));
    funtang1 = @(x) epsilon*(double(x<a));
elseif isempty(a)
    fun1 = @(x) epsilon*(max(x-b,0));
    funtang1 = @(x) epsilon*(double(x>b) );
end
fun = @(x) fun1(expand(x));
funtang = @(x) funtang1(expand(x));


dim =length(size(u));
n = size(u);

if ischarin('SEPSOLVER',varargin)
    solver = getcharin('SEPSOLVER',varargin);
else
    solver = SEPSOLVER(dim,varargin{:});
end
param = getparam(solver);


rm = u;
um = SEPMATRIX(dim);

erriter = zeros(1,param.maxiter);
errorder = zeros(1,param.maxorder);
errorderinf = zeros(1,param.maxorder);
valmin = [];
valmax = [];

switch param.errorindicator
    case 'residual'
        normref = normND(u);
    case 'reference'
        normref = normND(u);
end

for i=1:param.maxorder
    if param.display
        fprintf('order #%d \n',i)
    end
    
    U0 = seprand(n);
    U0 = normalizefuns(U0);
    alpha0=0;
    U = U0;
    
    for kkk=1:param.maxiter    %%%%  ITERATION  ALTERNEE
        for j=1:dim %%%% BOUCLE DIMENSION j
            S = prodsepdouble(U',rm,j);
            
            w = solve(N,S,...
                @(w) w+prodsepdouble(U',fun(expand(um+setfun(U,w,j))),j),...
                @(w) speye(n(j),n(j))+diag(prodsepdouble((U.*U)',funtang(expand(um+setfun(U,w,j))),j)));
            
            alpha = norm(w);
            U{1,j} = w/alpha;
            
        end %%%% END BOUCLE DIMENSION j
        
        erriter(kkk)=abs(alpha-alpha0)/(alpha+alpha0);
        if param.display
            fprintf('  iteration #%d - stagnation = %d\n',kkk,erriter(kkk))
        end
        
        if erriter(kkk)<param.itercrit
            break
        end
        alpha0=alpha;
        % U0 = U;
        
    end  %%%%  END ITERATION ALTERNEE
    
    um = um + alpha*U;
    
    if param.update>0 && i>1 && mod(i,param.updatestep)==0
        for kkk=1:param.update
            if param.display
                fprintf('update #%d, dim ',kkk)
            end
            um0=um;
            W = gathervectors(um);
            
            if param.ortho
                W = myorth(W);
            end
            
            for jj=1:length(param.updatedim)
                j = param.updatedim(jj);
                if param.display
                    fprintf('#%d ',j)
                end
                % W.alpha = ones(1,i);
                % WW = W'*W;
                % WW.F(:,j)= {ones(i,i)};
                %
                % S = prodsepdouble(W',u,j);
                %
                % M = timesblock(WW);
                % M = M +param.updateeps*eye(size(M,1));
                % Wj = (M\S')';
                % W.alpha = ones(1,i);
                % W.F{1,j}=Wj;
                
                W.alpha = ones(1,i);
                S = prodsepdouble(W',u,j);
                % S = assembleblock(Wu,j);
                WW = W'*W;
                WW.F(:,j)= {speye(n(j),n(j))};
                M = assembleblock(WW,j);
                
                Wj = solve(N,S(:),@(Wj) M*Wj + funblock(Wj,W,fun,j,n,i),...
                    @(Wj) M + funblocktang(Wj,W,funtang,j,n,i));
                Wj = reshape(Wj,n(j),i);
                W.F{1,j}=Wj;
                
                
                
            end
            um = splitvectors(W);
            um = normalizefuns(um);
            
            if kkk<param.update
                errup0= sqrt(min(um.alpha)^2/sum(um.alpha.^2));
                % stagn = abs(errup-errup0)/errup;
                stagn = full(norm(um-um0)/norm(um0));
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
        rm = u - expand(um);
    else
        rm = rm - expand(alpha*U);
    end
    
    
    valmin(i) = minND(expand(um));
    valmax(i) = maxND(expand(um));
    
    switch param.errorindicator
        case {'residual','reference'}
            errorder(i)=normND(rm)/normref;
            errorderinf(i)=normNDinf(rm)/normNDinf(expand(u));
        otherwise
            errorder(i)= sqrt(min(um.alpha)^2/sum(um.alpha.^2));
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
result.errorinf = errorderinf;
result.valmin = valmin ;
result.valmax = valmax;

function U = setfun(U,w,j)
U{1,j} = w;
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



function Dj = funblock(Wj,W,fun,j,n,m)
Wj = reshape(Wj,n(j),m);
W.F{1,j}=Wj;
W=splitvectors(W);
Dj = zeros(n(j),m);
F = fun(W);
for k=1:m
    Dj(:,k) = prodsepdouble(truncate(W,k)',F,j);
end
Dj = Dj(:);

return


function D2j = funblocktang(Wj,W,funtang,j,n,m)

Wj = reshape(Wj,n(j),m);
W.F{1,j}=Wj;
W=splitvectors(W);
D2j = sparse(n(j)*m,n(j)*m);
F = funtang(W);
for k=1:m
    for l=1:m
        repk = ((k-1)*n(j)+1:k*n(j));
        repl = ((l-1)*n(j)+1:l*n(j));
        Wk = truncate(W,k);Wk.alpha(:)=1;
        Wl = truncate(W,l);Wl.alpha(:)=1;
        D2j(repk,repl) = diag(prodsepdouble((Wk.*Wl)',F,j));
    end
end


return

function ab = prodsepdouble(a,b,j)
m = size(a{1,1},1);
if nargin==2
    j=0;
    ab = zeros(1,m);
else
    ab = zeros(size(b,j),m);
end


for ii=1:m
    ab(:,ii) = prodsepdoubleuni(a,b,j,ii) ;
end

return

function ab = prodsepdoubleuni(a,b,j,ii)

stemp = size(b);
dim = length(stemp);
ab = b;
for jjj=1:dim
    if jjj~=j
        dt = 1:dim;
        dt(jjj)=1;
        dt(1)=jjj;
        ab = permute(ab,dt);
        ab = a{1,jjj}(ii,:)*ab(:,:);
        stemp = stemp(dt);
        stemp(1)=size(ab,1);
        ab = reshape(ab,stemp);
        stemp = stemp(dt);
        ab = permute(ab,dt);
    end
end

if j~=0
    dt = 1:dim;
    dt(j)=[];
    dt=[j,dt];
    ab = permute(ab,dt);
end

ab = squeeze(ab);


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



