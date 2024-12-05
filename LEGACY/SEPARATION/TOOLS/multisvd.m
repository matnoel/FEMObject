function [u,result] = multisvd(b,varargin)
% function [u,result] =  multisvd(b,varargin)


dim =length(size(b));
n = size(b);

if ischarin('SEPSOLVER',varargin)
    solver = getcharin('SEPSOLVER',varargin);
else
    solver = SEPSOLVER(dim,varargin{:});
end
param = getparam(solver);


bu = b;
u = SEPMATRIX(dim);

erriter = zeros(1,param.maxiter);
if ~isfield(param,'maxorder')
    param.maxorder = max(n);
end
errorder = zeros(1,param.maxorder);

switch param.errorindicator
    case 'residual'
        normref = normND(b);
    case 'reference'
        if isempty(param.reference)
            param.reference = b;
        end
        normref = normND(param.reference);
end

valmin = [];
valmax = [];
errorderinf = zeros(1,param.maxorder);

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
            
            S = prodsepdouble(U',bu,j);
            
            
            U{1,j} = S;
            alpha = norm(U{1,j});
            U{1,j} = U{1,j}/alpha;
            
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
    
    if param.update>0 && i>1 && mod(i,param.updatestep)==0
        for kkk=1:param.update
            if param.display
                fprintf('update #%d, dim ',kkk)
            end
            u0=u;
            W = gathervectors(u);
            
            if param.ortho
                W = myorth(W);
            end
            
            for jj=1:length(param.updatedim)
                j = param.updatedim(jj);
                if param.display
                    fprintf('#%d ',j)
                end
                alpha = W.alpha;
                Wj = W.F{1,j};
                W.alpha = ones(1,i);
                WW = W'*W;
                WW.F(:,j)= {ones(i,i)};
                
                S = prodsepdouble(W',b,j);
                
                M = timesblock(WW);
                M = M +param.updateeps*eye(size(M,1));
                if param.updateeps>0 || condest(M)<1e15
                    Wj = (M\S')';
                    alpha = ones(1,i);
                else
                    fprintf('(NO)')
                end
                W.alpha = alpha;
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
        bu = b - expand(u);
    else
        bu = bu - expand(alpha*U);
    end
    
    
    if param.alphaupdate
        W = gathervectors(u);
        W.alpha = ones(1,i);
        WW = timesblock(W'*W);
        Wb = prodsepdouble(W',b);
        W.alpha = (WW\Wb')';
        u=splitvectors(W);
        bu = b - expand(u);
    end
    
    valmin(i) = minND(expand(u));
    valmax(i) = maxND(expand(u));
    
    
    switch param.errorindicator
        case 'residual'
            errorder(i)=normND(bu)/normref;
        case 'reference'
            if isa(param.reference,'SEPMATRIX')
                errorder(i)=normND(u-param.reference)/normref;
            else
                errorder(i)=normND(expand(u)-param.reference)/normref;
            end
            errorderinf(i)=normNDinf(expand(u)-expand(param.reference))/normNDinf(expand(param.reference));
            
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
result.errorinf = errorderinf;
result.valmin = valmin ;
result.valmax = valmax;

%
% function M = timesblock(A)
%
% i=1;
% Mi = A.F{i,1};
% for k=2:A.dim
% Mi = Mi.*A.F{i,k};
% end
% M = Mi*A.alpha(i);
%
% for i=2:A.m
% Mi = A.F{i,1};
% for k=2:A.dim
% Mi = Mi.*A.F{i,k};
% end
% M = M + Mi*A.alpha(i);
% end
%
% return
%
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


%
% function W = myorth(W,j)
%
%     if size(W.F,1)>1
%         error('regrouper les fonctions')
%     end
%     if nargin==1
%         j = 1:size(W.F,2);
%     end
%     for k=1:length(j)
%        W.F{j(k)} = mygram(W.F{j(k)});
%     end
%
% return
%
% function W = mygram(W)
%
% for i=1:size(W,2)
% w = W(:,i);
% w = w/sqrt(w'*w);
% w = w - W(:,1:i-1)*(W(:,1:i-1)'*w);
% n = sqrt(w'*w);
% if n<1e-15
%     warning('orthogonalisation critique')
% end
% w = w/n;
% W(:,i) = w;
% end
%
% return
%
%
%
