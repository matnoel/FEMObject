function [u,utilde] = solve_updatetucker(A,b,u,utilde,param)

if nargin<=4
   param.adjoint=0;
   param.seltuckerindices=0;
end

if param.seltuckerindices
    orthfun = 0;
else
    orthfun=1;
end


dim = getdim(A);
%siz = repmat(getm(u),1,dim);
W = getF(gathervectors(u));
W = cellfun(@full,W,'UniformOutput',0);

if orthfun==1
W = cellfun(@(C) mygram(C,1e-10),W,'UniformOutput',0);
end
siz = cellfun(@(C) size(C,2),W,'UniformOutput',0);
siz = [siz{:}];
if param.adjoint 
Wtilde = getF(gathervectors(utilde)); 
Wtilde = cellfun(@full,Wtilde,'UniformOutput',0);
if orthfun==1
Wtilde = cellfun(@(C) mygram((C)),Wtilde,'UniformOutput',0);
end
else
Wtilde = W;    
end
%Wtilde = cellfun(@ctranspose,Wilde,'UniformOutput',false);

Wb = getF(gathervectors(b));
Wb = cellfun(@full,Wb,'UniformOutput',0);
Wb = cellfun(@(A,B) mtimes(A',B),Wtilde,Wb,'UniformOutput',false);

WAW = cell(dim,1);
for k=1:dim
WAW{k} = A.F(:,k);    
end

for k=1:dim
WAW{k} = cellfun(@(C) mtimes(Wtilde{k}',C), WAW{k},'UniformOutput',false);%WAW{k} = cellmtimes(Wtilde{k},A.F{k});    
WAW{k} = cellfun(@(C) mtimes(C,W{k}), WAW{k},'UniformOutput',false);%WAW{k} = cellmtimes(WAW{k},W{k});
WAW{k} = cellfun(@(C) reshape(C,siz(k)^2,1),WAW{k},'UniformOutput',false);
WAW{k} = [WAW{k}{:}];
end
%balpha = tensor(tensordiag(b.alpha,b.dim));% pourrait être amelioré  
%balpha =sptensordiag(b.alpha,b.dim);% pourrait être amelioré  
%S = double(ttm(balpha,Wb));
%Aalpha = tensor(tensordiag(A.alpha,A.dim));% pourrait être amelioré    
%Aalpha = sptensordiag(A.alpha,A.dim);% pourrait être amelioré    
%M = double(ttm(Aalpha,WAW));
b.F = Wb(:)';
A.F = WAW(:)';
S = expand(splitvectors(b));
M = expand(splitvectors(A));

siz2 = [siz;siz];siz2=siz2(:)';
M = reshape(M,siz2);
perm = [1:2:2*dim,2:2:2*dim];
M = permute(M,perm);
M = reshape(M,prod(siz),prod(siz));
S = reshape(S,prod(siz),1);

if param.seltuckerindices==0
D = M\S;    
else
I = selmultiindices(siz,'diagonaldistance',param.seltuckerindiceslevel);
D = zeros(size(S));
D(I) = M(I,I)\S(I);
end
D = reshape(D,siz);

if dim==2
% %W.F{1}=W.F{1}*DD;
% %W.alpha = ones(1,i);
 [UU,DD,VV]=svd(D,0);
 W{1}=W{1}*UU;
 W{2}=W{2}*VV;
 u.F=W;
 u.m=size(W{1},2);
 u = splitvectors(u);
 u.alpha = diag(DD)';
 u.m = length(u.alpha);
% u = normalizefuns(u);
else
    % approximation de l'optimal
    % on peut augmenter à souhait la précision ou maxorder dans multisvd
    % attention, la précision est la précision sur la solution finale
D= double(D);    
[D,restemp] = multisvd(D,'maxorder',u.m,'tol',param.tol/100,'reference',D,'update',1,'display',0);
fprintf(' ... approximation of core tensor: rank = %d, error = %d\n',getm(D),restemp.error(getm(D)));
%D = cp_als(tensor(D),u.m);D = SEPMATRIX(D.U(:)',D.lambda');
D = gathervectors(D);
D = checkdim(D,dim);

u.F = cellfun(@mtimes,W,D.F,'UniformOutput',false);
u = splitvectors(u);
u.alpha = D.alpha;
u.m = length(u.alpha);
end

if param.adjoint==1
warning('pas prévu')
end


function D = checkdim(D,dim)

if size(D.F,2)<dim
    temp = cell(size(D.F,1),dim-size(D.F,2));
    temp(:)={ones(1,length(D.alpha))};
    D.F = [D.F,temp];
end

return

    %[UU,DD,VV]=svd(DDD);
% W.F{1}=W.F{1}*UU;
% W.F{2}=W.F{2}*VV;
% W.alpha = diag(DD)';
% Wtilde.F{1}=Wtilde.F{1}*DDD;




% 
% dim = getdim(A);
% i = getm(u);
% W = gathervectors(u);  
% W.alpha = 1;
% W = myorth(W,1:dim);    
% if param.adjoint 
% Wtilde = gathervectors(utilde); 
% Wtilde.alpha = 1;
% Wtilde = myorth(Wtilde,1:dim);
% else
% Wtilde = W;    
% end
% 
% WAW = Wtilde'*A*W;
% Wb =  Wtilde'*b;
% 
% 
% M = zeros(i*i,i*i);
% S = zeros(i,i);
% for kkk=1:size(Wb.F,1)
%     S = S + Wb.alpha(kkk)*Wb.F{kkk,1}*Wb.F{kkk,2}';
% end
% S = reshape(S,[i*i,1]);
% for kkk=1:size(WAW.F,1)
%     M = M + WAW.alpha(kkk)*WAW.F{kkk,1}(:)*WAW.F{kkk,2}(:)';
% end
% M = reshape(M,[i,i,i,i]);
% M = permute(M,[1,3,2,4]);
% M = reshape(M,[i*i,i*i]);
% W.alpha = ones(1,i);
% 
% DD = reshape(M\S,[i,i]);
% %W.F{1}=W.F{1}*DD;
% %W.alpha = ones(1,i);
% [UU,DD,VV]=svd(DD);
% W.F{1}=W.F{1}*UU;
% W.F{2}=W.F{2}*VV;
% W.alpha = diag(DD)';
% u=splitvectors(W);
% u = normalizefuns(u);
% 
% if param.adjoint 
% Wtilde.F{1}=Wtilde.F{1}*DD;
% utilde = splitvectors(Wtilde);
% utilde = normalizefuns(utilde);
% end
% 
% 
% 
% 
