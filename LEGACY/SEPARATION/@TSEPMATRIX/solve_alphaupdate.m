function [u,utilde] = solve_alphaupdate(A,b,u,utilde,param)
% function [u,utilde] = solve_alphaupdate(A,b,u,utilde,param)

if nargin==3 || isempty(utilde)
    utilde = u;
end
if nargin<=4
    param.adjoint=0;
end

dim = getdim(A);
%siz = getm(u);
W = getF(gathervectors(u));
W = cellfun(@full,W,'UniformOutput',0);
W = cellfun(@(C) mygram((C)),W,'UniformOutput',0);

siz = cellfun(@(C)size(C,2),W)';
if param.adjoint 
Wtilde = getF(gathervectors(utilde)); 
Wtilde = cellfun(@full,Wtilde,'UniformOutput',0);
Wtilde = cellfun(@(C) mygram((C)),Wtilde,'UniformOutput',0);
else
Wtilde = W;    
end
Wtilde = cellctranspose(Wtilde);


%Wtilde = cellfun(@ctranspose,Wtilde,'UniformOutput',0);

% dim = getdim(A);
% u = orth(u);
% siz = getm(u);
% W = getF(gathervectors(u));  
% if param.adjoint 
% utilde = orth(utilde);
% Wtilde = getF(gathervectors(utilde)); 
% else
% Wtilde = W;    
% end
% Wtilde = cellfun(@ctranspose,Wtilde,'UniformOutput',0);
% %Wtilde = cellctranspose(Wtilde);

Wb = getF(gathervectors(b));
Wb = cellfun(@full,Wb,'UniformOutput',0);
Wb = cellfun(@mtimes,Wtilde,Wb,'UniformOutput',false);

WAW = A.F;
for k=1:dim
WAW{k} = cellfun(@(C) mtimes(Wtilde{k},C), WAW{k},'UniformOutput',false);%WAW{k} = cellmtimes(Wtilde{k},A.F{k});    
WAW{k} = cellfun(@(C) mtimes(C,W{k}), WAW{k},'UniformOutput',false);%WAW{k} = cellmtimes(WAW{k},W{k});
WAW{k} = cellfun(@(C) reshape(C,siz(k)^2,1),WAW{k},'UniformOutput',false);
WAW{k} = [WAW{k}{:}];
end
S = double(ttm(b.alpha,Wb));
M = double(ttm(A.alpha,WAW));
siz2 = [siz;siz];siz2=siz2(:)';
M = reshape(M,siz2);
perm = [1:2:2*dim,2:2:2*dim];
M = permute(M,perm);
M = reshape(M,prod(siz),prod(siz));
S = reshape(S,prod(siz),1);

D = M\S;
D = reshape(D,siz);
u.alpha = tensor(D);
u.alpha = reshape(D,siz);
u.F = W;
u = splitvectors(u);
%u = hooi(u,getm(u));

if param.adjoint==1
warning('pas prévu')
utilde.F = Wtilde;
utilde = splitvectors(utilde);
end
    %[UU,DD,VV]=svd(DDD);
% W.F{1}=W.F{1}*UU;
% W.F{2}=W.F{2}*VV;
% W.alpha = diag(DD)';
% Wtilde.F{1}=Wtilde.F{1}*DDD;

