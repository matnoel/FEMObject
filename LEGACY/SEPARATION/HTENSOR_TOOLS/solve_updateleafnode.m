function u = solve_updateleafnode(A,b,u,t)

N=numel(u.U{t});
uAu=zeros(N);
Au=apply_mat_to_vec(A,u);

M = reduceinnerprodatnode(u,Au,t);
Mb = reduceinnerprodatnode(u,b,t);

ub=ttm(b.U{t},Mb{3},2);
ub=ub(:);

ns = size(u);
nm = size(A);
m = nm./ns;
k_A = size(A.U{t}, 2);
mu = find(u.dim2ind == t);

szuUt = size(u.U{t});

%AUt=A.U{t};
AUt = mat2cell(A.U{t},size(A.U{t},1),ones(1,k_A));
Ajj = cellfun(@(x) reshape(x,m(mu),ns(mu)),AUt,'uniformoutput',false);

for n=1:N
    U_right=zeros(szuUt);
    U_right(n)=1;

%      AuU = cell(1, k_A);
%      for jj=1:k_A
%          Ajj = reshape(AUt(:, jj), m(mu), ns(mu));
%          AuU{jj} = Ajj * U_right;
%      end
    AuU = cellfun(@(x) x*U_right,Ajj,'uniformoutput',false);

    AuUt = cell2mat(AuU);

%AuUt = ttm(AuUt,M{3},2); too slow!
    AuUt = AuUt*M{3}';
    uAu(n,:)=AuUt(:);
end

u.U{t}=reshape(uAu\ub,size(u.U{t}));

end
