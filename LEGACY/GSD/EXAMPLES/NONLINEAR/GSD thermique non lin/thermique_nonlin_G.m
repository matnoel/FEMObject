function Lam=thermique_nonlin_G(W,tol,S,PC,X,varargin)

thermique_nonlin_forms

s2 = size(W,2);
a1 = X{1}*(W'*a{S}(:,:)*W);
a2 = [];
a3 = zeros(s2,s2,s2,s2);
for k=1:s2
    a3(:,:,k,k) = W'*g{S}(W(:,k),W(:,k),:,:)*W;
    for j=k+1:s2
        a3(:,:,j,k) = W'*g{S}(W(:,k),W(:,j),:,:)*W;
        a3(:,:,k,j)=a3(:,:,j,k);  
    end    
end
a3 = X{2}*a3;


if nbsec==1
    b = X{3}*W'*l{S}(:);
else
    b = X{3}*W'*l{S}(:)+ X{4}*W'*l2{S}(:);    
end

Lam = GSDthermiquesolvestoreac(PC,tol,a1,a2,a3,b,varargin{:});


